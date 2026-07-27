"""
k–ℓ 湍流闭合 + Kantha & Clayson (2004) / Harcourt 风格的 Langmuir 效应。

动量通量采用 Lagrangian / Harcourt 形式：
```
τ = -⟨u'w'⟩ = νt ( ∂U/∂z + αs ∂Us/∂z )
```
`αs=1` 即按 Lagrangian 剪切 `∂(U+Us)/∂z` 混合（McWilliams et al. 2012; Reichl et al. 2016）。
这是复现“欧拉力几乎可忽略、UL≈Us”的关键：应力可由 Stokes 剪切承担，无需巨大欧拉剪切。

TKE 源项用总应力（而非仅欧拉剪切）写 Stokes 生产：
```
P   = τ · ∂U/∂z
P_S = E6 · τ · ∂Us/∂z
ε   = cε k^{3/2}/ℓ
νt  = cμ √k · ℓ
```
"""
struct KLStokesClosure{T<:AbstractFloat}
    κ::T
    cμ::T
    cε::T
    σk::T
    E6::T
    αs::T          # Stokes 剪切在动量通量中的权重；1 = 完全 Lagrangian
    ℓ_max::T
    channel::Bool
    k_min::T
    νt_min::T
end

function KLStokesClosure(; κ = 0.4, cμ = 0.09, cε = 0.166, σk = 1.0, E6 = 4.0,
                         αs = 1.0, ℓ_max = Inf, channel = true,
                         k_min = 1e-12, νt_min = 0.0)
    T = Float64
    return KLStokesClosure{T}(T(κ), T(cμ), T(cε), T(σk), T(E6), T(αs),
                              T(float(ℓ_max)), channel, T(k_min), T(νt_min))
end

function eddy_viscosity_faces!(νt_f::AbstractVector, k::AbstractVector, ℓ::AbstractVector,
                               clos::KLStokesClosure, grid::UniformColumnGrid)
    Nz = grid.Nz
    cμ = clos.cμ
    νt_f[1] = clos.νt_min
    νt_f[end] = clos.νt_min
    @inbounds for i in 2:Nz
        qℓ = 0.5 * (sqrt(max(k[i - 1], clos.k_min)) * ℓ[i - 1] +
                    sqrt(max(k[i], clos.k_min)) * ℓ[i])
        νt_f[i] = max(cμ * qℓ, clos.νt_min)
    end
    return νt_f
end

"""
欧拉剪切生产与 Stokes 生产（基于总应力 τ = νt (Uz + αs Usz)）。
"""
function tke_production(Uz::Real, Vz::Real, Usz::Real, Vsz::Real, νt::Real,
                        E6::Real, αs::Real = 1.0)
    τu = νt * (Uz + αs * Usz)
    τv = νt * (Vz + αs * Vsz)
    P = τu * Uz + τv * Vz
    PS = E6 * (τu * Usz + τv * Vsz)
    return P, PS
end

# 兼容旧签名
tke_production(Uz, Vz, Usz, Vsz, νt, E6) = tke_production(Uz, Vz, Usz, Vsz, νt, E6, 1.0)

function update_tke!(k::AbstractVector, U::AbstractVector, V::AbstractVector,
                     stokes::StokesDrift, νt_f::AbstractVector, ℓ::AbstractVector,
                     clos::KLStokesClosure, grid::UniformColumnGrid, forcing::Forcing,
                     dt::Real)
    Nz = grid.Nz
    dz = grid.dz
    ν = forcing.ν
    σk = clos.σk
    cε = clos.cε
    E6 = clos.E6
    αs = clos.αs
    k_min = clos.k_min

    νt_c = similar(k)
    @inbounds for i in 1:Nz
        νt_c[i] = 0.5 * (νt_f[i] + νt_f[i + 1])
    end

    Uz = similar(k)
    Vz = similar(k)
    @inbounds for i in 1:Nz
        if i == 1
            Uz[i] = (U[i + 1] - U[i]) / dz
            Vz[i] = (V[i + 1] - V[i]) / dz
        elseif i == Nz
            Uz[i] = (U[i] - U[i - 1]) / dz
            Vz[i] = (V[i] - V[i - 1]) / dz
        else
            Uz[i] = (U[i + 1] - U[i - 1]) / (2dz)
            Vz[i] = (V[i + 1] - V[i - 1]) / (2dz)
        end
    end

    rhs = similar(k)
    @inbounds for i in 1:Nz
        P, PS = tke_production(Uz[i], Vz[i], stokes.dusdz_c[i], stokes.dvsdz_c[i],
                               νt_c[i], E6, αs)
        rhs[i] = P + PS
    end

    Kf = similar(νt_f)
    @inbounds for i in eachindex(νt_f)
        Kf[i] = ν + νt_f[i] / σk
    end
    @inbounds for i in 1:Nz
        Km = Kf[i]
        Kp = Kf[i + 1]
        km = i == 1 ? k[i] : k[i - 1]
        kp = i == Nz ? k[i] : k[i + 1]
        flux_m = Km * (k[i] - km) / dz
        flux_p = Kp * (kp - k[i]) / dz
        rhs[i] += (flux_p - flux_m) / dz
    end

    @inbounds for i in 1:Nz
        √k = sqrt(max(k[i], k_min))
        denom = 1 + dt * cε * √k / max(ℓ[i], eps(typeof(ℓ[i])))
        k[i] = max((k[i] + dt * rhs[i]) / denom, k_min)
    end
    return k
end
