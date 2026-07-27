"""
K-profile 参数化 + Langmuir 增强（KPPLT）。

参考：
- Large et al. (1994) KPP
- McWilliams & Sullivan (2000); Li & Fox-Kemper (2017)

```
νt(σ) = h · w_s · G(σ),   G(σ)=σ(1-σ)^2
w_s = κ u★ √(1 + Cw / La_t²)     (中性 + Langmuir)
```

默认 `Cw` 按 Xuan & Shen (2025) Fig.2b 的 νt 峰值量级校准：
对 `La_t=0.3`，`max νt/(u★H)≈0.38` ⇒ `Cw≈3.6`（`G_max=4/27`）。
"""
struct KPPLTClosure{T<:AbstractFloat}
    κ::T
    Cw::T
    φm::T
    αs::T          # Lagrangian 应力权重
    νt_min::T
    use_langmuir::Bool
end

function KPPLTClosure(; κ = 0.4, Cw = 3.6, φm = 1.0, αs = 1.0,
                      νt_min = 0.0, use_langmuir = true)
    T = Float64
    return KPPLTClosure{T}(T(κ), T(Cw), T(φm), T(αs), T(νt_min), use_langmuir)
end

shape_G(σ::Real) = σ * (1 - σ)^2

function eddy_viscosity_kpp!(νt_f::AbstractVector, clos::KPPLTClosure,
                             grid::UniformColumnGrid, forcing::Forcing, La_t::Real;
                             h = nothing)
    H = grid.H
    hBL = isnothing(h) ? H : oftype(H, h)
    u★ = forcing.u★
    ws = clos.κ * u★ / clos.φm
    if clos.use_langmuir && La_t > 0
        ws *= sqrt(1 + clos.Cw / La_t^2)
    end

    @inbounds for i in eachindex(grid.zf)
        z = grid.zf[i]
        σ = clamp(-z / hBL, 0, 1)
        νt_f[i] = max(hBL * ws * shape_G(σ), clos.νt_min)
    end
    return νt_f
end
