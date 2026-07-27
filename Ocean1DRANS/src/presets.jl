"""
    xuan_shen_config(; Nz=64, H=1.0, La_t=0.3, Reτ=1000, k0H=3.5, ...)

Xuan & Shen (2025) 型无分层 Langmuir 通道：
- 表面常风应力，底部 stress-free
- 常压力梯度 `Fx = -u★²/H` 平衡风应力
- 无 Coriolis、无浮力
- 深水指数 Stokes：`Us = (u★/La_t²) exp(2 k0H z/H)`

默认返回**无量纲**配置（速度用 `u★`、长度用 `H` 无量纲化）。

关键字 `closure`：
- `:my25` / `:kc04`（默认）— Mellor–Yamada 2.5 + Kantha–Clayson (2004)；GOTM 对齐
- `:harcourt` — 同上，但动量用 Lagrangian 应力 (`αs=1`)
- `:les` — Fig.2b LES 数字化 νt；论文对照 / resolvent 基流
- `:kpplt` — 峰值校准的 KPPLT
- `:klstokes` — 简化代数 k–ℓ（仅作趋势对照，勿当成熟模型）
- 或传入具体闭合对象
"""
function xuan_shen_config(;
        Nz::Integer = 64,
        H::Real = 1.0,
        La_t::Real = 0.3,
        Reτ::Real = 1000.0,
        k0H::Real = 3.5,
        u★::Real = 1.0,
        nondimensional::Bool = true,
        closure = :my25,
        E6::Real = 4.0,
        αs::Real = 0.0,
        κ::Real = 0.4,
        cμ::Real = 0.09,
        cε::Real = 0.166,
        Cw::Real = 3.6,
        Sm0::Real = 0.39327,
    )
    T = Float64
    H = T(H)
    u★ = T(u★)
    ν = u★ * H / T(Reτ)
    Fx = -(u★^2) / H
    _ = nondimensional

    grid = UniformColumnGrid(Nz, H)
    forcing = Forcing(u★; f = 0.0, Fx = Fx, Fy = 0.0, ν = ν)
    stokes = monochromatic_stokes(grid; u★ = u★, La_t = La_t, k0H = k0H)
    boundary = BoundarySetup(bottom = :stress_free)
    initial = InitialState(U = 0.0, V = 0.0, k = nothing)

    clos = if closure isa Symbol
        if closure === :my25 || closure === :kc04
            MY25KC04Closure(; κ = κ, E6 = E6, αs = αs, Sm0 = Sm0)
        elseif closure === :harcourt
            HarcourtMomentumClosure(; κ = κ, E6 = E6, Sm0 = Sm0)
        elseif closure === :klstokes
            KLStokesClosure(; κ = κ, cμ = cμ, cε = cε, E6 = E6, αs = 1.0, channel = true)
        elseif closure === :kpplt
            KPPLTClosure(; κ = κ, Cw = Cw, αs = 1.0, use_langmuir = true)
        elseif closure === :les
            LESNutClosure(; La_t = La_t, αs = 1.0)
        else
            error("unknown closure symbol $closure (use :my25, :kc04, :harcourt, :klstokes, :kpplt, :les)")
        end
    else
        closure
    end

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
        closure = :my25,
        E6::Real = 4.0,
        αs::Real = 0.0,
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

    clos = if closure isa Symbol
        if closure === :my25 || closure === :kc04
            MY25KC04Closure(; E6 = E6, αs = αs)
        elseif closure === :harcourt
            HarcourtMomentumClosure(; E6 = E6)
        elseif closure === :klstokes
            KLStokesClosure(; E6 = E6, αs = 1.0, channel = false)
        else
            error("unknown closure symbol $closure")
        end
    elseif isnothing(closure)
        MY25KC04Closure(; E6 = E6, αs = αs)
    else
        closure
    end

    return ModelConfig{T,typeof(clos)}(
        grid, forcing, stokes, boundary, initial, clos,
        T(La_t), T(k0H), false,
    )
end
