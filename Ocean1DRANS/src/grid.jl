"""
均匀分层的一维水柱网格。

约定：
- `z` 向上为正，表面 `z=0`，底边界 `z=-H`
- 速度与 TKE 定义在层中心 `zc`
- 通量与涡粘定义在层界面 `zf`（含底/表共 `Nz+1` 个界面）
"""
struct UniformColumnGrid{T<:AbstractFloat}
    Nz::Int
    H::T
    dz::T
    zc::Vector{T}   # cell centers, length Nz (bottom → surface)
    zf::Vector{T}   # faces, length Nz+1 (bottom → surface)
end

function UniformColumnGrid(Nz::Integer, H::Real)
    Nz >= 4 || throw(ArgumentError("Nz must be >= 4"))
    H > 0 || throw(ArgumentError("H must be positive"))
    T = float(typeof(H))
    H = T(H)
    dz = H / Nz
    # faces from bottom to surface
    zf = collect(range(-H, 0; length = Nz + 1))
    zc = @. 0.5 * (zf[1:end-1] + zf[2:end])
    return UniformColumnGrid{T}(Int(Nz), H, dz, zc, zf)
end

cell_centers(g::UniformColumnGrid) = g.zc
interfaces(g::UniformColumnGrid) = g.zf
depths(g::UniformColumnGrid) = .-g.zc  # positive depth
