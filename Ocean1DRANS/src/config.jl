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

"""
    xuan_shen_config(; Nz=64, H=1.0, La_t=0.3, Reτ=1000, k0H=3.5, ...)

Xuan & Shen (2025) 型无分层 Langmuir 通道：
- 表面常风应力，底部 stress-free
- 常压力梯度 `Fx = -u★²/H` 平衡风应力
- 无 Coriolis、无浮力
- 深水指数 Stokes：`Us = (u★/La_t²) exp(2 k0H z/H)`

默认返回**无量纲**配置（速度用 `u★`、长度用 `H` 无量纲化），
此时 `u★=1, H=1, ν=1/Reτ`。设 `nondimensional=false` 可改用有量纲。
"""
function xuan_shen_config(;
        Nz::Integer = 64,
        H::Real = 1.0,
        La_t::Real = 0.3,
        Reτ::Real = 1000.0,
        k0H::Real = 3.5,
        u★::Real = 1.0,
        nondimensional::Bool = true,
        closure = nothing,
        E6::Real = 4.0,
        κ::Real = 0.4,
        cμ::Real = 0.09,
        cε::Real = 0.166,   # ≈ cμ^{3/4}/κ^{1/2} 量级；与 MOST 匹配时再调
    )
    T = Float64
    H = T(H)
    u★ = T(u★)
    # 无量纲模式下默认 H=1, u★=1，分子粘性由 Reτ = u★ H / ν 确定
    ν = u★ * H / T(Reτ)
    Fx = -(u★^2) / H
    _ = nondimensional  # 标记：调用方可据此解释输出量纲

    grid = UniformColumnGrid(Nz, H)
    forcing = Forcing(u★; f = 0.0, Fx = Fx, Fy = 0.0, ν = ν)
    stokes = monochromatic_stokes(grid; u★ = u★, La_t = La_t, k0H = k0H)
    boundary = BoundarySetup(bottom = :stress_free)
    initial = InitialState(U = 0.0, V = 0.0, k = nothing)

    clos = isnothing(closure) ?
        KLStokesClosure(; κ = κ, cμ = cμ, cε = cε, E6 = E6, channel = true) :
        closure

    return ModelConfig{T,typeof(clos)}(
        grid, forcing, stokes, boundary, initial, clos,
        T(La_t), T(k0H), nondimensional,
    )
end

"""
    mcwilliams1997_config(; Nz=50, H=90, u★=0.0061, La_t=0.3, ...)

McWilliams et al. (1997) 型开洋混合层：含 Coriolis 与 Stokes–Coriolis，
底边界二次拖曳，Stokes 波长默认 60 m。
"""
function mcwilliams1997_config(;
        Nz::Integer = 50,
        H::Real = 90.0,
        u★::Real = 0.0061,
        La_t::Real = 0.3,
        wavelength::Real = 60.0,
        f::Real = 1e-4,
        ν::Real = 1e-6,
        Cd::Real = 1.5e-3,
        closure = nothing,
        E6::Real = 4.0,
    )
    T = Float64
    H = T(H)
    u★ = T(u★)
    k0 = 2π / T(wavelength)
    k0H = k0 * H

    grid = UniformColumnGrid(Nz, H)
    forcing = Forcing(u★; f = T(f), Fx = 0.0, Fy = 0.0, ν = T(ν))
    stokes = monochromatic_stokes(grid; u★ = u★, La_t = La_t, k0H = k0H)
    boundary = BoundarySetup(bottom = :quadratic_drag, Cd = Cd)
    initial = InitialState(U = 0.0, V = 0.0)

    clos = isnothing(closure) ?
        KLStokesClosure(; E6 = E6, channel = false) :
        closure

    return ModelConfig{T,typeof(clos)}(
        grid, forcing, stokes, boundary, initial, clos,
        T(La_t), T(k0H), false,
    )
end
