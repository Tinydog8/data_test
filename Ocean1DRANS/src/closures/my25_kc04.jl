"""
Mellor–Yamada 2.5 / Kantha–Clayson (2004) / GOTM-style two-equation closure.

Prognostic variables: ``q²`` (twice TKE) and ``q²ℓ``.

KC04 Langmuir extension: Stokes production enters **both** equations, with
``E6 = 4`` in the length equation. Without that term, ``KM`` collapses
relative to LES (KC04 Fig.1).

References:
- Kantha & Clayson (2004), Ocean Modelling 6:101–124
- Mellor & Yamada (1982); GOTM Scientific Documentation
- Harcourt (2013, 2015) optional Lagrangian stress (`αs=1`)
"""

"""
    MY25KC04Closure(; E6=4.0, αs=0.0, ...)

Kantha–Clayson (2004) extension of Mellor–Yamada 2.5.

- `E6`: Stokes coefficient in ``q²ℓ`` (KC04 recommend 4)
- `αs`: momentum Stokes weight (0 = strict KC04; 1 = Harcourt Lagrangian)
- `Sm0`: equilibrium stability function ≈ 0.39327
"""
struct MY25KC04Closure{T<:AbstractFloat}
    B1::T
    Sq::T
    Sl::T
    E1::T
    E2::T
    E6::T
    κ::T
    Sm0::T
    αs::T
    q2_min::T
    ℓ_min::T
    νt_min::T
    νt_max::T
end

function MY25KC04Closure(; B1 = 16.6, Sq = 0.2, Sl = 0.2, E1 = 1.8, E2 = 1.33,
                         E6 = 4.0, κ = 0.4, Sm0 = 0.39327, αs = 0.0,
                         q2_min = 1e-10, ℓ_min = 1e-6, νt_min = 0.0, νt_max = Inf)
    T = Float64
    return MY25KC04Closure{T}(T(B1), T(Sq), T(Sl), T(E1), T(E2), T(E6), T(κ), T(Sm0),
                              T(αs), T(q2_min), T(ℓ_min), T(νt_min), T(νt_max))
end

"""Kantha–Clayson (2004) defaults (`E6=4`, `αs=0`)."""
KanthaClayson2004Closure(; kwargs...) = MY25KC04Closure(; kwargs...)

"""Deprecated: use [`Harcourt2015Closure`](@ref) for full SMC; this is KC04+`αs=1`."""
KC04LagrangianClosure(; kwargs...) = MY25KC04Closure(; αs = 1.0, kwargs...)
"""
Wall proximity length ``L_z = κ z_s z_b / (z_s + z_b)`` (MY / GOTM).
"""
function wall_length_Lz!(Lz::AbstractVector, grid::UniformColumnGrid, κ::Real)
    H = grid.H
    @inbounds for i in eachindex(Lz)
        zs = max(-grid.zc[i], eps(typeof(H)))          # distance to surface
        zb = max(H + grid.zc[i], eps(typeof(H)))       # distance to bottom
        Lz[i] = κ * zs * zb / (zs + zb)
    end
    return Lz
end

function initialize_my25!(q2::AbstractVector, q2l::AbstractVector,
                          grid::UniformColumnGrid, u★::Real, clos::MY25KC04Closure)
    Lz = similar(q2)
    wall_length_Lz!(Lz, grid, clos.κ)
    q2s = (clos.B1^(2 / 3)) * u★^2
    @inbounds for i in eachindex(q2)
        ℓ = max(0.1 * Lz[i], clos.ℓ_min)
        q2[i] = max(q2s * (Lz[i] / max(maximum(Lz), eps())), clos.q2_min)
        q2l[i] = q2[i] * ℓ
    end
    # Surface Dirichlet scale
    q2[end] = max(q2s, clos.q2_min)
    q2l[end] = q2[end] * max(clos.κ * (grid.dz / 2), clos.ℓ_min)
    return q2, q2l
end

"""
Eddy viscosity at cell centers: ``KM = q ℓ S_M``, then average to faces.
"""
function eddy_viscosity_my25!(νt_f::AbstractVector, νt_c::AbstractVector,
                              q2::AbstractVector, q2l::AbstractVector,
                              clos::MY25KC04Closure, grid::UniformColumnGrid)
    Nz = grid.Nz
    H = grid.H
    @inbounds for i in 1:Nz
        q2i = max(q2[i], clos.q2_min)
        ℓ = clamp(q2l[i] / q2i, clos.ℓ_min, 0.5 * H)
        νt_c[i] = clamp(sqrt(q2i) * ℓ * clos.Sm0, clos.νt_min, clos.νt_max)
    end
    νt_f[1] = clos.νt_min
    νt_f[end] = clos.νt_min
    @inbounds for i in 2:Nz
        νt_f[i] = 0.5 * (νt_c[i - 1] + νt_c[i])
    end
    return νt_f
end

function _shear_at_centers!(Uz::AbstractVector, Vz::AbstractVector,
                            U::AbstractVector, V::AbstractVector, dz::Real)
    Nz = length(U)
    @inbounds for i in 1:Nz
        if i == 1
            Uz[i] = (U[2] - U[1]) / dz
            Vz[i] = (V[2] - V[1]) / dz
        elseif i == Nz
            Uz[i] = (U[Nz] - U[Nz - 1]) / dz
            Vz[i] = (V[Nz] - V[Nz - 1]) / dz
        else
            Uz[i] = (U[i + 1] - U[i - 1]) / (2dz)
            Vz[i] = (V[i + 1] - V[i - 1]) / (2dz)
        end
    end
    return Uz, Vz
end

"""
Shear and Stokes production (KC04):
``P = KM S²``, ``P_s = KM (∂U/∂z ∂Us/∂z + ∂V/∂z ∂Vs/∂z)``.
"""
function my25_production(Uz::Real, Vz::Real, Usz::Real, Vsz::Real, KM::Real)
    P = KM * (Uz^2 + Vz^2)
    Ps = KM * (Uz * Usz + Vz * Vsz)
    return P, Ps
end

"""
Local algebraic equilibration of ``q²`` and ``q²ℓ`` (production–dissipation + wall),
then optional diffusion relaxation. Used inside the stress-balance steady solver.
"""
function equilibrate_my25!(q2::AbstractVector, q2l::AbstractVector,
                           U::AbstractVector, V::AbstractVector,
                           stokes::StokesDrift, νt_c::AbstractVector,
                           clos::MY25KC04Closure, grid::UniformColumnGrid, u★::Real;
                           underrelax::Real = 0.5, n_diff::Integer = 8)
    Nz = grid.Nz
    dz = grid.dz
    H = grid.H
    Lz = similar(q2)
    wall_length_Lz!(Lz, grid, clos.κ)
    Uz = similar(q2)
    Vz = similar(q2)
    _shear_at_centers!(Uz, Vz, U, V, dz)

    q2_new = similar(q2)
    q2l_new = similar(q2l)
    @inbounds for i in 1:Nz
        KM = max(νt_c[i], clos.νt_min)
        P, Ps = my25_production(Uz[i], Vz[i], stokes.dusdz_c[i], stokes.dvsdz_c[i], KM)
        # From q²: P + Ps = q³/(B1 ℓ)  ⇒  q³ = B1 ℓ (P+Ps)
        # From q²ℓ wall sink: E1 P + E6 Ps = (q³/B1) (1 + E2 (ℓ/Lz)²) / ℓ * ℓ
        #   = (q³/B1) (1 + E2 (ℓ/Lz)²)
        # Substitute q³/B1 = ℓ (P+Ps):
        #   E1 P + E6 Ps = ℓ (P+Ps) (1 + E2 (ℓ/Lz)²)
        # Solve for ℓ (scalar Newton / fixed-point).
        Prod = max(P + Ps, 0.0)
        Src_l = max(clos.E1 * P + clos.E6 * Ps, 0.0)
        if Prod < 1e-18 || Src_l < 1e-18
            q2_new[i] = clos.q2_min
            q2l_new[i] = clos.q2_min * clos.ℓ_min
            continue
        end
        Lzi = max(Lz[i], clos.ℓ_min)
        # Fixed point: ℓ = Src_l / (Prod * (1 + E2 (ℓ/Lz)²))
        ℓ = clamp(q2l[i] / max(q2[i], clos.q2_min), clos.ℓ_min, 0.5 * H)
        for _ in 1:12
            wall = 1 + clos.E2 * (ℓ / Lzi)^2
            ℓ_tgt = Src_l / (Prod * wall)
            ℓ = 0.5 * ℓ + 0.5 * clamp(ℓ_tgt, clos.ℓ_min, 0.5 * H)
        end
        # q³ = B1 ℓ Prod  ⇒  q² = (B1 ℓ Prod)^{2/3}
        q2_new[i] = max((clos.B1 * ℓ * Prod)^(2 / 3), clos.q2_min)
        q2l_new[i] = q2_new[i] * ℓ
    end

    # Surface Dirichlet (MY wall law)
    q2s = (clos.B1^(2 / 3)) * u★^2
    ℓs = max(clos.κ * (dz / 2), clos.ℓ_min)
    q2_new[end] = max(q2s, clos.q2_min)
    q2l_new[end] = q2_new[end] * ℓs

    @. q2 = underrelax * q2_new + (1 - underrelax) * q2
    @. q2l = underrelax * q2l_new + (1 - underrelax) * q2l
    q2[end] = max(q2s, clos.q2_min)
    q2l[end] = q2[end] * ℓs

    # Mild vertical diffusion of q², q²ℓ (Sq, Sl) for smoothness
    for _ in 1:n_diff
        _diffuse_scalar_explicit!(q2, νt_c, clos.Sq, grid, clos.q2_min)
        _diffuse_scalar_explicit!(q2l, νt_c, clos.Sl, grid, clos.q2_min * clos.ℓ_min)
        q2[end] = max(q2s, clos.q2_min)
        q2l[end] = q2[end] * ℓs
    end
    return q2, q2l
end

function _diffuse_scalar_explicit!(φ::AbstractVector, νt_c::AbstractVector, σ::Real,
                                   grid::UniformColumnGrid, φ_min::Real;
                                   cfl::Real = 0.2)
    Nz = grid.Nz
    dz = grid.dz
    φn = copy(φ)
    @inbounds for i in 2:Nz-1
        Km = 0.5 * (νt_c[i - 1] + νt_c[i]) * σ
        Kp = 0.5 * (νt_c[i] + νt_c[i + 1]) * σ
        Kmax = max(Km, Kp, eps(typeof(dz)))
        dt = cfl * dz^2 / Kmax
        flux_m = Km * (φ[i] - φ[i - 1]) / dz
        flux_p = Kp * (φ[i + 1] - φ[i]) / dz
        φn[i] = max(φ[i] + dt * (flux_p - flux_m) / dz, φ_min)
    end
    # bottom Neumann
    @inbounds begin
        Kp = 0.5 * (νt_c[1] + νt_c[2]) * σ
        dt = cfl * dz^2 / max(Kp, eps(typeof(dz)))
        flux_p = Kp * (φ[2] - φ[1]) / dz
        φn[1] = max(φ[1] + dt * flux_p / dz, φ_min)
    end
    φ .= φn
    return φ
end

"""
Time-advance ``q²`` and ``q²ℓ`` one step (for Coriolis / unsteady path).
"""
function advance_my25!(q2::AbstractVector, q2l::AbstractVector,
                       U::AbstractVector, V::AbstractVector,
                       stokes::StokesDrift, νt_c::AbstractVector, νt_f::AbstractVector,
                       clos::MY25KC04Closure, grid::UniformColumnGrid, u★::Real, dt::Real)
    Nz = grid.Nz
    dz = grid.dz
    H = grid.H
    Lz = similar(q2)
    wall_length_Lz!(Lz, grid, clos.κ)
    Uz = similar(q2)
    Vz = similar(q2)
    _shear_at_centers!(Uz, Vz, U, V, dz)

    rhs_q2 = zeros(eltype(q2), Nz)
    rhs_q2l = zeros(eltype(q2), Nz)
    @inbounds for i in 1:Nz
        q2i = max(q2[i], clos.q2_min)
        q = sqrt(q2i)
        ℓ = clamp(q2l[i] / q2i, clos.ℓ_min, 0.5 * H)
        KM = max(νt_c[i], clos.νt_min)
        P, Ps = my25_production(Uz[i], Vz[i], stokes.dusdz_c[i], stokes.dvsdz_c[i], KM)
        ε = q^3 / (clos.B1 * ℓ)
        rhs_q2[i] = 2 * (P + Ps) - 2 * ε
        wall = 1 + clos.E2 * (ℓ / max(Lz[i], clos.ℓ_min))^2
        rhs_q2l[i] = ℓ * (clos.E1 * P + clos.E6 * Ps) - (q^3 / clos.B1) * wall
    end

    # diffusion ∂z( Sq q ℓ ∂φ/∂z ) ≈ ∂z( Sq KM/Sm0 * ... ); use Sq * qℓ ≈ Kq
    Kq_f = similar(νt_f)
    @inbounds for i in eachindex(νt_f)
        Kq_f[i] = clos.Sq / max(clos.Sm0, eps()) * νt_f[i]
    end
    @inbounds for i in 1:Nz
        Km = Kq_f[i]
        Kp = Kq_f[i + 1]
        if i == 1
            flux_m_q2 = 0.0
            flux_m_q2l = 0.0
        else
            flux_m_q2 = Km * (q2[i] - q2[i - 1]) / dz
            flux_m_q2l = Km * (q2l[i] - q2l[i - 1]) / dz
        end
        if i == Nz
            flux_p_q2 = 0.0
            flux_p_q2l = 0.0
        else
            flux_p_q2 = Kp * (q2[i + 1] - q2[i]) / dz
            flux_p_q2l = Kp * (q2l[i + 1] - q2l[i]) / dz
        end
        rhs_q2[i] += (flux_p_q2 - flux_m_q2) / dz
        rhs_q2l[i] += (flux_p_q2l - flux_m_q2l) / dz
    end

    @inbounds for i in 1:Nz
        q2[i] = max(q2[i] + dt * rhs_q2[i], clos.q2_min)
        q2l[i] = max(q2l[i] + dt * rhs_q2l[i], clos.q2_min * clos.ℓ_min)
    end
    # Surface Dirichlet
    q2s = (clos.B1^(2 / 3)) * u★^2
    ℓs = max(clos.κ * (dz / 2), clos.ℓ_min)
    q2[end] = max(q2s, clos.q2_min)
    q2l[end] = q2[end] * ℓs
    return q2, q2l
end

function sync_tke_from_q2!(k::AbstractVector, ℓ::AbstractVector,
                           q2::AbstractVector, q2l::AbstractVector,
                           clos::MY25KC04Closure, H::Real)
    @inbounds for i in eachindex(k)
        q2i = max(q2[i], clos.q2_min)
        k[i] = 0.5 * q2i          # TKE = q²/2
        ℓ[i] = clamp(q2l[i] / q2i, clos.ℓ_min, 0.5 * H)
    end
    return k, ℓ
end
