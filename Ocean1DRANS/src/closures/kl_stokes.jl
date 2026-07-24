"""
k–ℓ 湍流闭合 + Kantha & Clayson (2004) / GOTM 风格的 Langmuir（Stokes）生产项。

TKE 方程（稳态/时间推进）：
```
∂k/∂t = P + E6 * P_S + ∂z[(ν + νt/σk) ∂z k] - ε
```
其中
```
P   = νt (Uz² + Vz²)                 # 欧拉剪切生产
P_S = νt (Uz Usz + Vz Vsz)           # Stokes 剪切生产
ε   = cε k^{3/2} / ℓ
νt  = cμ √k · ℓ
```

`E6 > 1`（KC04 建议约 4）用于增强 Langmuir 混合对涡粘的贡献。
"""
struct KLStokesClosure{T<:AbstractFloat}
    κ::T
    cμ::T
    cε::T
    σk::T
    E6::T          # Stokes 生产放大系数 (Kantha & Clayson 2004)
    ℓ_max::T
    channel::Bool  # true: 通道型双壁混合长度
    k_min::T
    νt_min::T
end

function KLStokesClosure(; κ = 0.4, cμ = 0.09, cε = 0.166, σk = 1.0, E6 = 4.0,
                         ℓ_max = Inf, channel = true, k_min = 1e-12, νt_min = 0.0)
    T = Float64
    return KLStokesClosure{T}(T(κ), T(cμ), T(cε), T(σk), T(E6), T(float(ℓ_max)),
                              channel, T(k_min), T(νt_min))
end

function eddy_viscosity_faces!(νt_f::AbstractVector, k::AbstractVector, ℓ::AbstractVector,
                               clos::KLStokesClosure, grid::UniformColumnGrid)
    # interpolate √k * ℓ to faces
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

function tke_production(Uz::Real, Vz::Real, Usz::Real, Vsz::Real, νt::Real, E6::Real)
    P = νt * (Uz^2 + Vz^2)
    PS = νt * (Uz * Usz + Vz * Vsz)
    return P, E6 * PS
end

"""
一步更新 TKE（向前欧拉 + 隐式耗散线性化）。
`νt_c` 为层中心涡粘（由界面平均得到）。
"""
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
    k_min = clos.k_min

    # cell-center viscosity for production
    νt_c = similar(k)
    @inbounds for i in 1:Nz
        νt_c[i] = 0.5 * (νt_f[i] + νt_f[i+1])
    end

    # shear at cell centers (central difference using neighboring cells;
    # one-sided near boundaries via face reconstruction)
    Uz = similar(k)
    Vz = similar(k)
    @inbounds for i in 1:Nz
        if i == 1
            Uz[i] = (U[i+1] - U[i]) / dz
            Vz[i] = (V[i+1] - V[i]) / dz
        elseif i == Nz
            Uz[i] = (U[i] - U[i-1]) / dz
            Vz[i] = (V[i] - V[i-1]) / dz
        else
            Uz[i] = (U[i+1] - U[i-1]) / (2dz)
            Vz[i] = (V[i+1] - V[i-1]) / (2dz)
        end
    end

    # RHS without dissipation (explicit) + diffusion
    rhs = similar(k)
    @inbounds for i in 1:Nz
        P, PS = tke_production(Uz[i], Vz[i], stokes.dusdz_c[i], stokes.dvsdz_c[i],
                               νt_c[i], E6)
        rhs[i] = P + PS
    end

    # vertical diffusion of k: ∂z[(ν + νt/σk) ∂z k]
    # face diffusivity
    Kf = similar(νt_f)
    @inbounds for i in eachindex(νt_f)
        Kf[i] = ν + νt_f[i] / σk
    end
    # zero-flux BCs for k at top/bottom → ghost via Neumann
    @inbounds for i in 1:Nz
        Km = Kf[i]
        Kp = Kf[i+1]
        km = i == 1 ? k[i] : k[i-1]
        kp = i == Nz ? k[i] : k[i+1]
        flux_m = Km * (k[i] - km) / dz
        flux_p = Kp * (kp - k[i]) / dz
        rhs[i] += (flux_p - flux_m) / dz
    end

    # implicit dissipation: k_new = (k + dt*rhs) / (1 + dt * cε √k / ℓ)
    @inbounds for i in 1:Nz
        √k = sqrt(max(k[i], k_min))
        denom = 1 + dt * cε * √k / max(ℓ[i], eps(typeof(ℓ[i])))
        k[i] = max((k[i] + dt * rhs[i]) / denom, k_min)
    end
    return k
end
