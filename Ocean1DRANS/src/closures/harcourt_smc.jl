"""
Harcourt (2015) second-moment closure for Langmuir turbulence.

Ported from GOTM `cmue_d_h15.F90` / `surface_proximity_function` /
`production` (STOKESFLUX) / `lengthscaleeq_kc04` (H15 ``S_l``).

Momentum flux (ARSM):
```
⟨u'w'⟩ = -qℓ ( S_M ∂U/∂z + S_S ∂Us/∂z )
```
i.e. ``KM = S_M q ℓ``, ``K_M^S = S_S q ℓ``.

TKE production (consistent with fluxes):
```
P  = KM |∂U/∂z|² + K_M^S (∂U/∂z·∂Us/∂z)
P_S = KM (∂U/∂z·∂Us/∂z) + K_M^S |∂Us/∂z|²
```

Surface proximity SPF = 1 - f_z^S suppresses near-surface Stokes
contributions (H15 fix for H13 defects). Default ``E6 = 6`` (H15).
"""

struct Harcourt2015Closure{T<:AbstractFloat}
    A1::T
    A2::T
    B1::T
    B2::T
    C1::T
    C2::T
    C3::T
    Sq::T
    E1::T
    E2::T
    E6::T
    κ::T
    ClS::T          # surface-proximity constant (H15 / GOTM: 0.25)
    Ghmin::T
    Ghoff::T
    Gvoff::T
    Sxmax::T
    q2_min::T
    ℓ_min::T
    νt_min::T
    νt_max::T
end

function Harcourt2015Closure(; A1 = 0.92, A2 = 0.74, B1 = 16.6, B2 = 10.1,
                             C1 = 0.08, C2 = 0.7, C3 = 0.2,
                             Sq = 0.2, E1 = 1.8, E2 = 1.33, E6 = 6.0,
                             κ = 0.4, ClS = 0.25,
                             Ghmin = -0.28, Ghoff = 0.003, Gvoff = 0.006,
                             Sxmax = 2.12,
                             q2_min = 1e-10, ℓ_min = 1e-6,
                             νt_min = 0.0, νt_max = Inf)
    T = Float64
    return Harcourt2015Closure{T}(
        T(A1), T(A2), T(B1), T(B2), T(C1), T(C2), T(C3),
        T(Sq), T(E1), T(E2), T(E6), T(κ), T(ClS),
        T(Ghmin), T(Ghoff), T(Gvoff), T(Sxmax),
        T(q2_min), T(ℓ_min), T(νt_min), T(νt_max))
end

# Keep old name as alias to the full SMC (no longer the αs=1 KC04 hack)
HarcourtMomentumClosure(; kwargs...) = Harcourt2015Closure(; kwargs...)

"""H15 ``S_l`` (GOTM lengthscaleeq_kc04 Eq.37): ``√(0.04+(0.41 S_H)²)``."""
harcourt_Sl(Sh::Real) = sqrt(0.04 + (0.41 * Sh)^2)

function initialize_harcourt!(q2::AbstractVector, q2l::AbstractVector,
                              grid::UniformColumnGrid, u★::Real, clos::Harcourt2015Closure)
    Lz = similar(q2)
    wall_length_Lz!(Lz, grid, clos.κ)
    q2s = (clos.B1^(2 / 3)) * u★^2
    @inbounds for i in eachindex(q2)
        ℓ = max(0.1 * Lz[i], clos.ℓ_min)
        q2[i] = max(q2s * (Lz[i] / max(maximum(Lz), eps())), clos.q2_min)
        q2l[i] = q2[i] * ℓ
    end
    q2[end] = max(q2s, clos.q2_min)
    q2l[end] = q2[end] * max(clos.κ * (grid.dz / 2), clos.ℓ_min)
    return q2, q2l
end

"""
Surface proximity SPF = 1 - f_z^S, ``f_z^S = 1 + tanh(C_l^S z / ℓ^S)``,
``ℓ^S`` = dissipation length weighted by positive Stokes production.
``z`` upward, surface at 0 (GOTM convention).
"""
function surface_proximity_SPF!(SPF::AbstractVector, ℓ::AbstractVector,
                                PS::AbstractVector, grid::UniformColumnGrid,
                                clos::Harcourt2015Closure)
    Nz = grid.Nz
    dz = grid.dz
    lS_num = 0.0
    lS_den = 0.0
    @inbounds for i in 1:Nz
        Psp = max(PS[i], 0.0)
        lS_num += ℓ[i] * Psp * dz
        lS_den += Psp * dz
    end
    lS = lS_num / max(lS_den, 1e-14)
    lS = max(lS, clos.ℓ_min)
    @inbounds for i in 1:Nz
        # zz = zc (surface 0, bottom -H)
        fzS = 1 + tanh(clos.ClS * grid.zc[i] / lS)
        SPF[i] = 1 - fzS
        SPF[i] = clamp(SPF[i], 0.0, 1.0)
    end
    return SPF
end

"""
Harcourt (2015) / GOTM `cmue_d_h15` stability functions for one cell.
Returns ``(S_M, S_S, S_H)`` in Mellor–Yamada (qℓ) convention.
"""
function harcourt_stability(Gh::Real, Gm::Real, Gv::Real, Gs::Real, SPF::Real,
                            clos::Harcourt2015Closure; length_lim::Bool = false)
    A1, A2, B1, B2 = clos.A1, clos.A2, clos.B1, clos.B2
    C1, C2, C3 = clos.C1, clos.C2, clos.C3
    small = 1e-8

    Shn0 = A2 * (1 - 6 * A1 / B1)
    Shnh = -9 * A1 * A2 * (A2 * (1 - 6 * A1 / B1))
    Shns = 9 * A1 * A2 * (1 - 6 * A1 / B1) * (2 * A1 + A2)
    Shnv = 9 * A1 * A2 * (A2 * (1 - 6 * A1 / B1 - 3 * C1) -
                          2 * A1 * (1 - 6 * A1 / B1 + 3 * C1))
    Shdah = -9 * A1 * A2
    Shdav = -36 * A1 * A1
    Shdbh = -3 * A2 * (6 * A1 + B2 * (1 - C3))
    Shdv = -9 * A2 * A2 * (1 - C2)
    Shdvh = -162 * A1 * A1 * A2 * (2 * A1 + (2 - C2) * A2)
    Shdvv = 324 * A1 * A1 * A2 * A2 * (1 - C2)
    Ssn0 = A1 * (1 - 6 * A1 / B1)
    Ssdh = -9 * A1 * A2
    Ssdv = -9 * A1 * A1
    Smn0 = A1 * (1 - 6 * A1 / B1 - 3 * C1)
    SmnhSh = 9 * A1 * (2 * A1 + A2 * (1 - C2))
    SmnsSs = 27 * A1 * A1
    Smdh = -9 * A1 * A2
    Smdv = -36 * A1 * A1

    # Stable-side Galperin-like limit inside ARSM when length_lim=false
    if !length_lim
        tmp1 = 1.0
        tmp2 = clos.Ghmin / min(clos.Ghmin, Gh)
        tmp1 = min(tmp1, tmp2)
        if tmp1 < 1
            Gh *= tmp1
            Gv *= tmp1
            Gs *= tmp1
        end
    end

    # Unstable-side limiter (GOTM cmue_d_h15)
    tmp0 = 2.0
    if Gv > 0
        tmp1 = (Shdah + Shdbh) * Gh + (Shdav + Shdv) * Gv
        tmp1 += (Shdah * clos.Ghoff + Shdav * clos.Gvoff) * (Shdbh * Gh) +
                (Shdvh * clos.Ghoff + Shdvv * clos.Gvoff) * Gv
        tmp1 += (Shdah * Gh + Shdav * Gv) * (Shdbh * clos.Ghoff) +
                (Shdvh * Gh + Shdvv * Gv) * clos.Gvoff
        tmp2 = (Shdah * Gh + Shdav * Gv) * (Shdbh * Gh) +
               (Shdvh * Gh + Shdvv * Gv) * Gv
        tmp4 = 1 + (Shdah + Shdbh) * clos.Ghoff + (Shdav + Shdv) * clos.Gvoff +
               (Shdah * clos.Ghoff + Shdav * clos.Gvoff) * (Shdbh * clos.Ghoff) +
               (Shdvh * clos.Ghoff + Shdvv * clos.Gvoff) * clos.Gvoff
        tmp3 = tmp1 * tmp1 - 4 * tmp2 * tmp4
        if tmp3 >= 0 && tmp2 < 0
            tmp3 = (-tmp1 + sqrt(tmp3)) / (2 * tmp2)
        elseif tmp3 >= 0 && tmp2 > 0
            tmp3 = (-tmp1 - sqrt(tmp3)) / (2 * tmp2)
        else
            tmp3 = 2.0
        end
        if 0 < tmp3 < 1
            tmp0 = min(tmp0, tmp3)
        end
    end

    # Apply SPF after first limiter
    Gv *= SPF
    Gs *= SPF^2

    if Gh > 0
        tmp1 = 2 * (Shdah + Shdbh) * Gh + (Shdav + Shdv) * Gv
        tmp2 = (2 * Shdah * Gh + Shdav * Gv) * (2 * Shdbh * Gh) +
               (2 * Shdvh * Gh + Shdvv * Gv) * Gv
        tmp4 = 1.0
        tmp3 = tmp1 * tmp1 - 4 * tmp2 * tmp4
        if tmp3 >= 0 && tmp2 < 0
            tmp3 = (-tmp1 + sqrt(tmp3)) / (2 * tmp2)
        elseif tmp3 >= 0 && tmp2 > 0
            tmp3 = (-tmp1 - sqrt(tmp3)) / (2 * tmp2)
        else
            tmp3 = 2.0
        end
        if 0 < tmp3 < 1
            tmp0 = min(tmp0, tmp3)
        end
    end

    if 0 < tmp0 < 1
        Gh *= tmp0
        Gm *= tmp0
        Gv *= tmp0
        Gs *= tmp0
    end

    # Sh
    num = Shn0 + Shnh * Gh + Shns * Gs + Shnv * Gv
    if num < 0
        Sh = small
    else
        den = (1 + Shdah * Gh + Shdav * Gv) * (1 + Shdbh * Gh) +
              (Shdv + Shdvh * Gh + Shdvv * Gv) * Gv
        Sh = den <= 0 ? clos.Sxmax : clamp(num / den, small, clos.Sxmax)
    end

    # Ss
    den = 1 + Ssdh * Gh + Ssdv * Gv
    Ss = den < 0 ? clos.Sxmax : clamp(Ssn0 / den, small, clos.Sxmax)

    # Sm
    num = Smn0 + SmnhSh * Gh * Sh + SmnsSs * Gs * Ss
    if abs(num) < small
        Gh += copysign(small, num == 0 ? 1.0 : num)
        Gv += copysign(small, num == 0 ? 1.0 : num)
        num = Smn0 + SmnhSh * Gh * Sh + SmnsSs * Gs * Ss
    end
    if num < 0
        Sm = small
    else
        den = 1 + Smdh * Gh + Smdv * Gv
        Sm = den <= 0 ? clos.Sxmax : clamp(num / den, small, clos.Sxmax)
    end

    Ss *= SPF
    return Sm, Ss, Sh
end

"""
Nondimensional forcings in MY variables, then H15 stability.
``G_M = (ℓ/q)² S²``, related to GOTM ``α`` by ``G = (4/B1²) α``.
"""
function harcourt_stability_from_state(Uz::Real, Vz::Real, Usz::Real, Vsz::Real,
                                       q2::Real, ℓ::Real, SPF::Real,
                                       clos::Harcourt2015Closure)
    q2 = max(q2, clos.q2_min)
    q = sqrt(q2)
    ℓ = max(ℓ, clos.ℓ_min)
    # ε = q³/(B1 ℓ); GOTM τ = k/ε with k=q²/2 ⇒ τ² S² = (B1²/4) G_M
    # Work directly in G_* :
    lq2 = (ℓ / q)^2
    S2 = Uz^2 + Vz^2
    CS = Uz * Usz + Vz * Vsz
    SS = Usz^2 + Vsz^2
    Gm = lq2 * S2
    Gv = lq2 * CS
    Gs = lq2 * SS
    Gh = 0.0   # unstratified default; stratified path can pass N² later
    return harcourt_stability(Gh, Gm, Gv, Gs, SPF, clos)
end

function harcourt_diffusivities!(KM::AbstractVector, KMS::AbstractVector,
                                 Sm::AbstractVector, Ss::AbstractVector, Sh::AbstractVector,
                                 q2::AbstractVector, q2l::AbstractVector,
                                 U::AbstractVector, V::AbstractVector,
                                 stokes::StokesDrift, SPF::AbstractVector,
                                 clos::Harcourt2015Closure, grid::UniformColumnGrid)
    Nz = grid.Nz
    dz = grid.dz
    H = grid.H
    Uz = similar(q2)
    Vz = similar(q2)
    _shear_at_centers!(Uz, Vz, U, V, dz)
    @inbounds for i in 1:Nz
        q2i = max(q2[i], clos.q2_min)
        ℓ = clamp(q2l[i] / q2i, clos.ℓ_min, 0.5 * H)
        q = sqrt(q2i)
        Sm[i], Ss[i], Sh[i] = harcourt_stability_from_state(
            Uz[i], Vz[i], stokes.dusdz_c[i], stokes.dvsdz_c[i],
            q2i, ℓ, SPF[i], clos)
        KM[i] = clamp(Sm[i] * q * ℓ, clos.νt_min, clos.νt_max)
        KMS[i] = clamp(Ss[i] * q * ℓ, clos.νt_min, clos.νt_max)
    end
    return KM, KMS
end

"""Consistent Harcourt production from fluxes."""
function harcourt_production(Uz::Real, Vz::Real, Usz::Real, Vsz::Real,
                             KM::Real, KMS::Real)
    S2 = Uz^2 + Vz^2
    CS = Uz * Usz + Vz * Vsz
    SS = Usz^2 + Vsz^2
    P = KM * S2 + KMS * CS
    Ps = KM * CS + KMS * SS
    return P, Ps
end

function faces_from_centers!(ν_f::AbstractVector, ν_c::AbstractVector, ν_min::Real)
    Nz = length(ν_c)
    ν_f[1] = ν_min
    ν_f[end] = ν_min
    @inbounds for i in 2:Nz
        ν_f[i] = 0.5 * (ν_c[i - 1] + ν_c[i])
    end
    return ν_f
end

"""
Algebraic equilibration of ``q²``, ``q²ℓ`` with H15 production and ``S_l(S_H)``.
Inner iterations close the ARSM–TKE coupling at fixed velocity.
"""
function equilibrate_harcourt!(q2::AbstractVector, q2l::AbstractVector,
                               U::AbstractVector, V::AbstractVector,
                               stokes::StokesDrift,
                               KM::AbstractVector, KMS::AbstractVector, Sh::AbstractVector,
                               clos::Harcourt2015Closure, grid::UniformColumnGrid, u★::Real;
                               underrelax::Real = 0.5, n_inner::Integer = 8,
                               n_diff::Integer = 2)
    Nz = grid.Nz
    dz = grid.dz
    H = grid.H
    Lz = similar(q2)
    wall_length_Lz!(Lz, grid, clos.κ)
    Uz = similar(q2)
    Vz = similar(q2)
    _shear_at_centers!(Uz, Vz, U, V, dz)

    SPF = ones(eltype(q2), Nz)
    Sm = similar(q2)
    Ss = similar(q2)
    PS = similar(q2)
    ℓ_tmp = similar(q2)

    q2s = (clos.B1^(2 / 3)) * u★^2
    ℓs = max(clos.κ * (dz / 2), clos.ℓ_min)

    for _inner in 1:n_inner
        @inbounds for i in 1:Nz
            ℓ_tmp[i] = clamp(q2l[i] / max(q2[i], clos.q2_min), clos.ℓ_min, 0.5 * H)
            _, Ps = harcourt_production(Uz[i], Vz[i], stokes.dusdz_c[i], stokes.dvsdz_c[i],
                                        KM[i], KMS[i])
            PS[i] = Ps
        end
        surface_proximity_SPF!(SPF, ℓ_tmp, PS, grid, clos)
        harcourt_diffusivities!(KM, KMS, Sm, Ss, Sh, q2, q2l, U, V, stokes, SPF, clos, grid)

        @inbounds for i in 1:Nz
            P, Ps = harcourt_production(Uz[i], Vz[i], stokes.dusdz_c[i], stokes.dvsdz_c[i],
                                        KM[i], KMS[i])
            Prod = max(P + Ps, 0.0)
            Src_l = max(clos.E1 * P + clos.E6 * Ps, 0.0)
            if Prod < 1e-18 || Src_l < 1e-18
                q2_tgt = clos.q2_min
                q2l_tgt = clos.q2_min * clos.ℓ_min
            else
                Lzi = max(Lz[i], clos.ℓ_min)
                ℓ = ℓ_tmp[i]
                for _ in 1:8
                    wall = 1 + clos.E2 * (ℓ / Lzi)^2
                    ℓ_tgt = Src_l / (Prod * wall)
                    ℓ = 0.5 * ℓ + 0.5 * clamp(ℓ_tgt, clos.ℓ_min, 0.5 * H)
                end
                q2_tgt = max((clos.B1 * ℓ * Prod)^(2 / 3), clos.q2_min)
                q2l_tgt = q2_tgt * ℓ
            end
            q2[i] = underrelax * q2_tgt + (1 - underrelax) * q2[i]
            q2l[i] = underrelax * q2l_tgt + (1 - underrelax) * q2l[i]
        end
        q2[end] = max(q2s, clos.q2_min)
        q2l[end] = q2[end] * ℓs
    end

    # Final diffusivity consistent with relaxed state
    @inbounds for i in 1:Nz
        ℓ_tmp[i] = clamp(q2l[i] / max(q2[i], clos.q2_min), clos.ℓ_min, 0.5 * H)
        _, Ps = harcourt_production(Uz[i], Vz[i], stokes.dusdz_c[i], stokes.dvsdz_c[i],
                                    KM[i], KMS[i])
        PS[i] = Ps
    end
    surface_proximity_SPF!(SPF, ℓ_tmp, PS, grid, clos)
    harcourt_diffusivities!(KM, KMS, Sm, Ss, Sh, q2, q2l, U, V, stokes, SPF, clos, grid)

    for _ in 1:n_diff
        Kq = similar(q2)
        @inbounds for i in 1:Nz
            q = sqrt(max(q2[i], clos.q2_min))
            ℓ = clamp(q2l[i] / max(q2[i], clos.q2_min), clos.ℓ_min, 0.5 * H)
            Sl = harcourt_Sl(Sh[i])
            Kq[i] = max(clos.Sq, 0.5 * Sl) * q * ℓ
        end
        _diffuse_scalar_explicit!(q2, Kq, 1.0, grid, clos.q2_min)
        _diffuse_scalar_explicit!(q2l, Kq, 1.0, grid, clos.q2_min * clos.ℓ_min)
        q2[end] = max(q2s, clos.q2_min)
        q2l[end] = q2[end] * ℓs
    end
    harcourt_diffusivities!(KM, KMS, Sm, Ss, Sh, q2, q2l, U, V, stokes, SPF, clos, grid)
    return q2, q2l, SPF
end

function sync_tke_from_q2_h!(k::AbstractVector, ℓ::AbstractVector,
                             q2::AbstractVector, q2l::AbstractVector,
                             clos::Harcourt2015Closure, H::Real)
    @inbounds for i in eachindex(k)
        q2i = max(q2[i], clos.q2_min)
        k[i] = 0.5 * q2i
        ℓ[i] = clamp(q2l[i] / q2i, clos.ℓ_min, 0.5 * H)
    end
    return k, ℓ
end

"""
Time-advance ``q²``, ``q²ℓ`` one step (Coriolis / unsteady path).
"""
function advance_harcourt!(q2::AbstractVector, q2l::AbstractVector,
                           U::AbstractVector, V::AbstractVector,
                           stokes::StokesDrift,
                           KM::AbstractVector, KMS::AbstractVector, Sh::AbstractVector,
                           SPF::AbstractVector,
                           clos::Harcourt2015Closure, grid::UniformColumnGrid,
                           u★::Real, dt::Real)
    Nz = grid.Nz
    dz = grid.dz
    H = grid.H
    Lz = similar(q2)
    wall_length_Lz!(Lz, grid, clos.κ)
    Uz = similar(q2)
    Vz = similar(q2)
    _shear_at_centers!(Uz, Vz, U, V, dz)

    PS = similar(q2)
    @inbounds for i in 1:Nz
        _, Ps = harcourt_production(Uz[i], Vz[i], stokes.dusdz_c[i], stokes.dvsdz_c[i],
                                    KM[i], KMS[i])
        PS[i] = Ps
    end
    ℓ_tmp = similar(q2)
    @inbounds for i in 1:Nz
        ℓ_tmp[i] = clamp(q2l[i] / max(q2[i], clos.q2_min), clos.ℓ_min, 0.5 * H)
    end
    surface_proximity_SPF!(SPF, ℓ_tmp, PS, grid, clos)
    Sm = similar(q2)
    Ss = similar(q2)
    harcourt_diffusivities!(KM, KMS, Sm, Ss, Sh, q2, q2l, U, V, stokes, SPF, clos, grid)

    rhs_q2 = zeros(eltype(q2), Nz)
    rhs_q2l = zeros(eltype(q2), Nz)
    @inbounds for i in 1:Nz
        q2i = max(q2[i], clos.q2_min)
        q = sqrt(q2i)
        ℓ = clamp(q2l[i] / q2i, clos.ℓ_min, 0.5 * H)
        P, Ps = harcourt_production(Uz[i], Vz[i], stokes.dusdz_c[i], stokes.dvsdz_c[i],
                                    KM[i], KMS[i])
        ε = q^3 / (clos.B1 * ℓ)
        rhs_q2[i] = 2 * (P + Ps) - 2 * ε
        wall = 1 + clos.E2 * (ℓ / max(Lz[i], clos.ℓ_min))^2
        rhs_q2l[i] = ℓ * (clos.E1 * P + clos.E6 * Ps) - (q^3 / clos.B1) * wall
    end

    @inbounds for i in 1:Nz
        q = sqrt(max(q2[i], clos.q2_min))
        ℓ = clamp(q2l[i] / max(q2[i], clos.q2_min), clos.ℓ_min, 0.5 * H)
        Sl = harcourt_Sl(Sh[i])
        Kqi = max(clos.Sq, 0.5 * Sl) * q * ℓ
        # store temporarily in Sm as diffusivity for compact diffusion
        Sm[i] = Kqi
    end
    Kq_f = zeros(eltype(q2), Nz + 1)
    faces_from_centers!(Kq_f, Sm, 0.0)
    @inbounds for i in 1:Nz
        Km = Kq_f[i]
        Kp = Kq_f[i + 1]
        flux_m = i == 1 ? 0.0 : Km * (q2[i] - q2[i - 1]) / dz
        flux_p = i == Nz ? 0.0 : Kp * (q2[i + 1] - q2[i]) / dz
        rhs_q2[i] += (flux_p - flux_m) / dz
        flux_m = i == 1 ? 0.0 : Km * (q2l[i] - q2l[i - 1]) / dz
        flux_p = i == Nz ? 0.0 : Kp * (q2l[i + 1] - q2l[i]) / dz
        rhs_q2l[i] += (flux_p - flux_m) / dz
    end

    @inbounds for i in 1:Nz
        q2[i] = max(q2[i] + dt * rhs_q2[i], clos.q2_min)
        q2l[i] = max(q2l[i] + dt * rhs_q2l[i], clos.q2_min * clos.ℓ_min)
    end
    q2s = (clos.B1^(2 / 3)) * u★^2
    ℓs = max(clos.κ * (dz / 2), clos.ℓ_min)
    q2[end] = max(q2s, clos.q2_min)
    q2l[end] = q2[end] * ℓs
    return q2, q2l
end
