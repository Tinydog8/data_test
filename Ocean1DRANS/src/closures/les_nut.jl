"""
直接使用 Xuan & Shen (2025) Fig.2(b) LES 数字化涡粘廓线的诊断闭合。

这是与论文 resolvent 分析最一致的 `νt` 来源：欧拉力可忽略时，
背景流取 `UL ≈ Us`，涡粘取 LES 剖面。
"""
struct LESNutClosure{T<:AbstractFloat}
    La_t::T
    αs::T
    y::Vector{T}       # y/H ∈ [0, -1]
    nu_t::Vector{T}    # νt / (u★ H)
end

function _default_les_csv()
    # src/closures -> package root/data
    return joinpath(@__DIR__, "..", "..", "data", "les_eddy_viscosity_fig2b.csv")
end

function LESNutClosure(; La_t::Real = 0.3, αs::Real = 1.0,
                       csv_path::AbstractString = _default_les_csv())
    T = Float64
    y, n02, n03 = load_les_nut_csv(csv_path)
    nu = isapprox(La_t, 0.2; atol = 1e-6) ? n02 :
         isapprox(La_t, 0.3; atol = 1e-6) ? n03 :
         error("LESNutClosure only tabulated for La_t = 0.2 or 0.3 (got $La_t)")
    return LESNutClosure{T}(T(La_t), T(αs), T.(y), T.(nu))
end

function eddy_viscosity_les!(νt_f::AbstractVector, clos::LESNutClosure,
                             grid::UniformColumnGrid, forcing::Forcing)
    u★H = forcing.u★ * grid.H
    @inbounds for i in eachindex(grid.zf)
        yH = grid.zf[i] / grid.H
        νt_f[i] = max(_interp_linear(clos.y, clos.nu_t, yH) * u★H, zero(u★H))
    end
    # 表面/底面保持小值（与应力边界协调）
    νt_f[1] = min(νt_f[1], νt_f[2])
    νt_f[end] = min(νt_f[end], νt_f[end - 1])
    return νt_f
end
