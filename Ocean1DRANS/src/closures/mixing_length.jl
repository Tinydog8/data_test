"""
混合长度诊断。

- 通道型（`channel=true`）：上下边界均限制，适用于 Xuan–Shen 型有限深通道
  ``ℓ = κ ζ_w √(1 - ζ_w/H_eff)``，其中 `ζ_w` 为到最近壁面的距离。
- 海洋混合层（`channel=false`）：Blackadar 型，仅由表面距离与 `ℓ_max` 限制
  ``ℓ^{-1} = (κ d)^{-1} + ℓ_max^{-1}``。
"""
function mixing_length!(ℓ::AbstractVector, grid::UniformColumnGrid;
                        κ::Real = 0.4,
                        ℓ_max::Real = Inf,
                        channel::Bool = false,
                        d_wall_min::Real = 1e-4)
    H = grid.H
    κ = oftype(H, κ)
    ℓ_max = oftype(H, ℓ_max)
    dmin = oftype(H, d_wall_min) * H

    if channel
        for i in eachindex(grid.zc)
            z = grid.zc[i]
            # distance to nearest boundary (surface z=0 or bottom z=-H)
            d = min(-z, z + H)
            d = max(d, dmin)
            # parabolic-channel-like length (Nezu & Rodi / open-channel form)
            ℓ[i] = κ * d * sqrt(max(one(H) - d / H, zero(H)))
            ℓ[i] = max(ℓ[i], κ * dmin)
        end
    else
        for i in eachindex(grid.zc)
            d = max(-grid.zc[i], dmin)
            if isfinite(ℓ_max)
                ℓ[i] = 1 / (1 / (κ * d) + 1 / ℓ_max)
            else
                ℓ[i] = κ * d
            end
        end
    end
    return ℓ
end

function mixing_length(grid::UniformColumnGrid; kwargs...)
    ℓ = similar(grid.zc)
    return mixing_length!(ℓ, grid; kwargs...)
end
