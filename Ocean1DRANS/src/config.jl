"""
外强迫参数。

- `u★`       : 表面摩擦速度 (m/s)
- `τx, τy`   : 表面动量通量 (= -ρ u★² 方向)；若为 `nothing` 则默认沿 -x：`τx=-u★², τy=0`
- `f`        : Coriolis 参数 (1/s)；通道型算例可取 0
- `Fx, Fy`   : 常体力/压力梯度 (m/s²)；Xuan–Shen 通道取 `Fx = -u★²/H`
- `ν`        : 分子粘性 (m²/s)
- `ρ0`       : 参考密度（仅文档用途）
"""
Base.@kwdef struct Forcing{T<:AbstractFloat}
    u★::T
    τx::T
    τy::T
    f::T = zero(T)
    Fx::T = zero(T)
    Fy::T = zero(T)
    ν::T = T(1e-6)
    ρ0::T = T(1025)
end

function Forcing(u★::Real; τx = nothing, τy = nothing, f = 0.0, Fx = 0.0, Fy = 0.0,
                 ν = 1e-6, ρ0 = 1025.0)
    T = float(typeof(float(u★)))
    u★ = T(u★)
    # 方程形式 ∂U/∂t = ∂z[(ν+νt)∂z U] + ...，表面通量 (ν+νt)∂zU = τx。
    # τx = +u★² 向流体注入动量；配合 Fx = -u★²/H 时底应力为 0。
    τx_val = isnothing(τx) ? u★^2 : T(τx)
    τy_val = isnothing(τy) ? zero(T) : T(τy)
    return Forcing{T}(u★, τx_val, τy_val, T(f), T(Fx), T(Fy), T(ν), T(ρ0))
end

"""
边界条件设置。

表面：通量型风应力（由 `Forcing` 给出）。
底部：
- `:stress_free` — 零切应力（Xuan–Shen）
- `:no_slip`     — U=V=0
- `:quadratic_drag` — 二次底拖曳，系数 `Cd`
"""
Base.@kwdef struct BoundarySetup
    bottom::Symbol = :stress_free
    Cd::Float64 = 1.5e-3
end

"""
初始条件。可传入函数 `z -> value` 或常数。
"""
Base.@kwdef struct InitialState
    U = 0.0
    V = 0.0
    k = nothing   # nothing → 由 u★ 估计
end

"""
完整模型配置。
"""
struct ModelConfig{T,C}
    grid::UniformColumnGrid{T}
    forcing::Forcing{T}
    stokes::StokesDrift{T}
    boundary::BoundarySetup
    initial::InitialState
    closure::C
    La_t::T
    k0H::T
    nondimensional::Bool
end
