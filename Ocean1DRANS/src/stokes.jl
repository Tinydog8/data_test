"""
Stokes 漂移廓线及其垂直剪切。

存储在层中心（与速度同位）与界面上，便于动量与湍流闭合同时使用。
"""
struct StokesDrift{T<:AbstractFloat}
    us_c::Vector{T}   # cell centers
    vs_c::Vector{T}
    us_f::Vector{T}   # faces
    vs_f::Vector{T}
    dusdz_c::Vector{T}
    dvsdz_c::Vector{T}
    dusdz_f::Vector{T}
    dvsdz_f::Vector{T}
end

"""
    exponential_stokes(grid; Us0, k0, angle=0)

深水单色波指数型 Stokes 漂移：
`uˢ(z) = Us0 * exp(2 k0 z)`（`z≤0`）。

Xuan & Shen (2025) / McWilliams et al. (1997) 常用形式。
`angle` 为相对 +x 轴的波向角（弧度）。
"""
function exponential_stokes(grid::UniformColumnGrid; Us0::Real, k0::Real, angle::Real = 0.0)
    T = typeof(grid.H)
    Us0 = T(Us0)
    k0 = T(k0)
    ca, sa = cos(T(angle)), sin(T(angle))

    us_c = similar(grid.zc)
    vs_c = similar(grid.zc)
    dusdz_c = similar(grid.zc)
    dvsdz_c = similar(grid.zc)
    for i in eachindex(grid.zc)
        amp = Us0 * exp(2k0 * grid.zc[i])
        damp = 2k0 * amp
        us_c[i] = ca * amp
        vs_c[i] = sa * amp
        dusdz_c[i] = ca * damp
        dvsdz_c[i] = sa * damp
    end

    us_f = similar(grid.zf)
    vs_f = similar(grid.zf)
    dusdz_f = similar(grid.zf)
    dvsdz_f = similar(grid.zf)
    for i in eachindex(grid.zf)
        amp = Us0 * exp(2k0 * grid.zf[i])
        damp = 2k0 * amp
        us_f[i] = ca * amp
        vs_f[i] = sa * amp
        dusdz_f[i] = ca * damp
        dvsdz_f[i] = sa * damp
    end

    return StokesDrift{T}(us_c, vs_c, us_f, vs_f, dusdz_c, dvsdz_c, dusdz_f, dvsdz_f)
end

"""
    monochromatic_stokes(grid; u★, La_t, k0H, angle=0)

由摩擦速度 `u★`、湍流 Langmuir 数 `La_t = √(u★/Us0)` 与无量纲波数 `k0H`
构造指数 Stokes 漂移。`Us0 = u★ / La_t^2`。
"""
function monochromatic_stokes(grid::UniformColumnGrid;
                              u★::Real, La_t::Real, k0H::Real = 3.5, angle::Real = 0.0)
    La_t > 0 || throw(ArgumentError("La_t must be positive"))
    Us0 = u★ / La_t^2
    k0 = k0H / grid.H
    return exponential_stokes(grid; Us0 = Us0, k0 = k0, angle = angle)
end
