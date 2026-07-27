"""
    xuan_shen_config(; Nz=64, H=1.0, La_t=0.3, Reτ=1000, k0H=3.5, ...)

Xuan & Shen (2025) 型无分层 Langmuir 通道：
- 表面常风应力，底部 stress-free
- 常压力梯度 `Fx = -u★²/H` 平衡风应力
- 无 Coriolis、无浮力
- 深水指数 Stokes：`Us = (u★/La_t²) exp(2 k0H z/H)`

默认返回**无量纲**配置（速度用 `u★`、长度用 `H` 无量纲化）。

关键字 `closure`：
- `:harcourt` / `:h15`（推荐）— Harcourt (2015) 完整 SMC（GOTM `cmue_d_h15`）
- `:my25` / `:kc04` — Mellor–Yamada 2.5 + Kantha–Clayson (2004)
- `:les` — Fig.2b LES 数字化 νt；论文对照 / resolvent 基流
- `:kpplt` — 峰值校准的 KPPLT
- `:klstokes` — 简化代数 k–ℓ（仅作趋势对照）
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
        closure = :harcourt,
        E6 = nothing,
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
        if closure === :harcourt || closure === :h15
            E6h = isnothing(E6) ? 6.0 : Float64(E6)
            Harcourt2015Closure(; κ = κ, E6 = E6h)
        elseif closure === :my25 || closure === :kc04
            E6m = isnothing(E6) ? 4.0 : Float64(E6)
            MY25KC04Closure(; κ = κ, E6 = E6m, αs = αs, Sm0 = Sm0)
        elseif closure === :kc04_lag
            E6m = isnothing(E6) ? 4.0 : Float64(E6)
            KC04LagrangianClosure(; κ = κ, E6 = E6m, Sm0 = Sm0)
        elseif closure === :klstokes
            E6m = isnothing(E6) ? 4.0 : Float64(E6)
            KLStokesClosure(; κ = κ, cμ = cμ, cε = cε, E6 = E6m, αs = 1.0, channel = true)
        elseif closure === :kpplt
            KPPLTClosure(; κ = κ, Cw = Cw, αs = 1.0, use_langmuir = true)
        elseif closure === :les
            LESNutClosure(; La_t = La_t, αs = 1.0)
        else
            error("unknown closure symbol $closure (use :harcourt, :my25, :kc04, :les, :kpplt, :klstokes)")
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
    mcwilliams1997_config(; Nz=64, H=33, u★=0.0061, La_t=0.3, ...)

McWilliams et al. (1997) / KC04 Fig.1 型开洋混合层：
- Coriolis（默认 45°，`f = 1.031×10⁻⁴`）与 Stokes–Coriolis
- 默认 `H = zi = 33 m`（McWilliams 反转层深度；域即混合层）
- 底边界默认 stress-free（混合层底近似，不当固壁）
- 单色波 λ=60 m，`La_t=0.3`，`u★=0.0061`
- MY25：`wall_mode=:surface`；`E6` 默认 **7.2**（Kantha et al. 2010 对 KC04
  原文 E6=4 笔误的更正；与 Fig.1 粗红线量级一致）

稳态求解走 Ekman 固定点 + 预后 q²/q²ℓ（非时间推进惯性振荡）。
"""
function mcwilliams1997_config(;
        Nz::Integer = 64,
        H::Real = 33.0,
        u★::Real = 0.0061,
        La_t::Real = 0.3,
        wavelength::Real = 60.0,
        f::Real = 1.031e-4,
        ν::Real = 1e-6,
        Cd::Real = 1.5e-3,
        bottom::Symbol = :stress_free,
        closure = :harcourt,
        E6 = nothing,
        αs::Real = 0.0,
        stokes_production::Bool = true,
    )
    T = Float64
    H = T(H)
    u★ = T(u★)
    k0 = 2π / T(wavelength)
    k0H = k0 * H

    grid = UniformColumnGrid(Nz, H)
    forcing = Forcing(u★; f = T(f), Fx = 0.0, Fy = 0.0, ν = T(ν))
    stokes = monochromatic_stokes(grid; u★ = u★, La_t = La_t, k0H = k0H)
    boundary = BoundarySetup(bottom = bottom, Cd = Cd)
    initial = InitialState(U = 0.0, V = 0.0)

    clos = if closure isa Symbol
        if closure === :harcourt || closure === :h15
            E6h = isnothing(E6) ? 6.0 : Float64(E6)
            Harcourt2015Closure(; E6 = E6h)
        elseif closure === :my25 || closure === :kc04
            # Kantha et al. (2010): KC04 printed E6=4 incorrectly; physical value ≈ 7.2
            E6m = isnothing(E6) ? 7.2 : Float64(E6)
            MY25KC04Closure(; E6 = E6m, αs = αs, stokes_production = stokes_production,
                            wall_mode = :ml, ℓ_max_frac = 0.6)
        elseif closure === :klstokes
            E6m = isnothing(E6) ? 4.0 : Float64(E6)
            KLStokesClosure(; E6 = E6m, αs = 1.0, channel = false)
        else
            error("unknown closure symbol $closure")
        end
    elseif isnothing(closure)
        Harcourt2015Closure()
    else
        closure
    end

    return ModelConfig{T,typeof(clos)}(
        grid, forcing, stokes, boundary, initial, clos,
        T(La_t), T(k0H), false,
    )
end
