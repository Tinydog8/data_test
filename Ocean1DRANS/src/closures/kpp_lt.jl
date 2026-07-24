"""
K-profile 参数化 + Langmuir 增强（KPPLT）。

参考：
- Large et al. (1994) KPP
- McWilliams & Sullivan (2000); Li & Fox-Kemper (2017) Langmuir 增强速度尺度

在无分层风生边界层中：
```
νt(σ) = h · w_s · G(σ),   σ = -z/h ∈ [0,1]
G(σ) = σ (1-σ)^2
w_s = κ u★ / φm ,   φm ≈ 1  (中性)
```
Langmuir 增强（LF17 简化）：
```
w_s ← w_s · √(1 + C_w / La_t²)
```
对通道型算例取 `h = H`；对开洋混合层可由 Ri 判据估计，此处默认 `h = H`
（无分层时混合至全深）。
"""
struct KPPLTClosure{T<:AbstractFloat}
    κ::T
    Cw::T       # Li & Fox-Kemper (2017) 量级常数
    φm::T       # 中性 MOST 稳定度函数
    νt_min::T
    use_langmuir::Bool
end

function KPPLTClosure(; κ = 0.4, Cw = 0.15, φm = 1.0, νt_min = 0.0, use_langmuir = true)
    T = Float64
    return KPPLTClosure{T}(T(κ), T(Cw), T(φm), T(νt_min), use_langmuir)
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
        # 表面/底界面涡粘置小值，避免与应力 BC 冲突时可保持光滑
        νt_f[i] = max(hBL * ws * shape_G(σ), clos.νt_min)
    end
    # 保证表面附近有足够混合：σ→0 时 G→0，对风应力驱动是合理的（近壁粘性次层由 ν 承担）
    return νt_f
end
