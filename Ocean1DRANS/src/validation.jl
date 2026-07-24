"""
物理一致性检查：应力平衡、Stokes 解析式、Langmuir 趋势、
TKE 生产–耗散平衡，以及与 LES 数字化涡粘剖面的对比。
"""

struct CheckResult
    name::String
    passed::Bool
    detail::String
    metric::Float64
end

function _pass(name, ok::Bool, detail::AbstractString; metric::Real = NaN)
    return CheckResult(name, ok, string(detail), Float64(metric))
end

# ---------------------------------------------------------------------------
# 1. 应力平衡（Xuan–Shen 通道）
# ---------------------------------------------------------------------------
function check_stress_balance(sol::SteadySolution; atol = 1e-10)
    cfg = sol.config
    diag = diagnostic_stress_balance(sol)
    τb = diag.τ_f[1]
    τt = diag.τ_f[end]
    ok_bottom = abs(τb) < atol
    ok_top = abs(τt - cfg.forcing.τx) < atol
    ok_pg = abs(cfg.forcing.τx + cfg.forcing.Fx * cfg.grid.H) < atol

    # 由 νe 与离散剪切重建界面应力，应贴近解析 τ(z)
    g = cfg.grid
    U = sol.state.U
    νe = diag.νe
    err = 0.0
    count = 0
    @inbounds for i in 2:g.Nz
        # 界面 i 介于 cell i-1 与 i
        Uz = (U[i] - U[i - 1]) / g.dz
        τ_num = νe[i] * Uz
        τ_ana = diag.τ_f[i]
        err += abs(τ_num - τ_ana)
        count += 1
    end
    mean_err = err / max(count, 1)
    ok_recon = mean_err < 1e-6 * max(abs(cfg.forcing.τx), 1.0)

    ok = ok_bottom && ok_top && ok_pg && ok_recon
    detail = @sprintf("τ_bottom=%.3e τ_top=%.3e  recon_MAE=%.3e", τb, τt, mean_err)
    return _pass("stress balance (τ, PG, reconstruction)", ok, detail; metric = mean_err)
end

# ---------------------------------------------------------------------------
# 2. Stokes 漂移解析式与 La_t 定义
# ---------------------------------------------------------------------------
function check_stokes_analytics(cfg::ModelConfig; rtol = 1e-12)
    g = cfg.grid
    Us0 = cfg.forcing.u★ / cfg.La_t^2
    k0 = cfg.k0H / g.H
    max_rel = 0.0
    @inbounds for i in eachindex(g.zc)
        ana = Us0 * exp(2k0 * g.zc[i])
        rel = abs(cfg.stokes.us_c[i] - ana) / max(abs(ana), eps())
        max_rel = max(max_rel, rel)
    end
    La_chk = sqrt(cfg.forcing.u★ / max(cfg.stokes.us_f[end], eps()))
    ok_la = abs(La_chk - cfg.La_t) / cfg.La_t < 1e-10
    ok = max_rel < rtol && ok_la
    detail = @sprintf("max|Us-Us_ana|/Us=%.3e  La_t(from Us0)=%.6f (target %.6f)",
                      max_rel, La_chk, cfg.La_t)
    return _pass("Stokes exponential profile & La_t definition", ok, detail; metric = max_rel)
end

# ---------------------------------------------------------------------------
# 3. 速度剪切符号与单调性（风应力驱动：∂U/∂z ≥ 0）
# ---------------------------------------------------------------------------
function check_velocity_structure(sol::SteadySolution)
    g = sol.config.grid
    U = sol.state.U
    # 表面速度应大于底部（φ[1]=0 规范）
    ok_surface = U[end] > U[1]
    # 绝大部分界面剪切应为正
    npos = 0
    ntot = g.Nz - 1
    @inbounds for i in 1:ntot
        if (U[i + 1] - U[i]) / g.dz >= -1e-12
            npos += 1
        end
    end
    frac = npos / ntot
    ok = ok_surface && frac > 0.95
    detail = @sprintf("U_surface=%.4f  U_bottom=%.4f  positive-shear fraction=%.3f",
                      U[end], U[1], frac)
    return _pass("Eulerian shear structure (wind-driven)", ok, detail; metric = frac)
end

# ---------------------------------------------------------------------------
# 4. TKE 局部生产–耗散平衡（KLStokes）
# ---------------------------------------------------------------------------
function check_tke_production_dissipation(sol::SteadySolution; rtol = 0.05)
    cfg = sol.config
    clos = cfg.closure
    clos isa KLStokesClosure || return _pass(
        "TKE production–dissipation balance", true,
        "skipped (closure is not KLStokes)"; metric = 0.0)

    g = cfg.grid
    U, V, k, ℓ = sol.state.U, sol.state.V, sol.state.k, sol.state.ℓ
    νt_c = sol.state.νt_c
    dz = g.dz
    max_rel = 0.0
    mean_rel = 0.0
    n = 0
    # 排除近壁 15% 厚度（混合长度迅速变化，差分剪切误差更大）
    @inbounds for i in 1:g.Nz
        σ = -g.zc[i] / g.H
        (σ < 0.15 || σ > 0.85) && continue
        if i == 1 || i == g.Nz
            continue
        end
        Uz = (U[i + 1] - U[i - 1]) / (2dz)
        Vz = (V[i + 1] - V[i - 1]) / (2dz)
        P, PS = Ocean1DRANS.tke_production(Uz, Vz, cfg.stokes.dusdz_c[i],
                                           cfg.stokes.dvsdz_c[i], νt_c[i], clos.E6)
        Prod = P + PS
        ε = clos.cε * (max(k[i], clos.k_min)^(3 / 2)) / max(ℓ[i], eps())
        if Prod > 1e-12
            rel = abs(Prod - ε) / Prod
            max_rel = max(max_rel, rel)
            mean_rel += rel
            n += 1
        end
    end
    mean_rel = n == 0 ? 0.0 : mean_rel / n
    ok = n > 0 && max_rel < rtol
    detail = @sprintf("max|P+PS-ε|/Prod=%.3e  mean=%.3e  (n=%d interior)", max_rel, mean_rel, n)
    return _pass("TKE production–dissipation balance", ok, detail; metric = max_rel)
end

# ---------------------------------------------------------------------------
# 5. Langmuir 增强混合：E6 与 La_t 趋势
# ---------------------------------------------------------------------------
function check_langmuir_mixing_trends(; Nz = 48)
    # E6: Stokes 生产开启应提高 νt
    sol_on = run_to_steady(xuan_shen_config(; Nz = Nz, La_t = 0.2, E6 = 4.0);
                           tol = 1e-6, max_steps = 500, verbose = false)
    sol_off = run_to_steady(xuan_shen_config(; Nz = Nz, La_t = 0.2, E6 = 0.0);
                            tol = 1e-6, max_steps = 500, verbose = false)
    ν_on = maximum(sol_on.state.νt_c)
    ν_off = maximum(sol_off.state.νt_c)
    ok_E6 = sol_on.converged && sol_off.converged && ν_on > ν_off

    # 更小 La_t（更强波）→ 更大 νt
    sol_02 = run_to_steady(xuan_shen_config(; Nz = Nz, La_t = 0.2, E6 = 4.0);
                           tol = 1e-6, max_steps = 500, verbose = false)
    sol_03 = run_to_steady(xuan_shen_config(; Nz = Nz, La_t = 0.3, E6 = 4.0);
                           tol = 1e-6, max_steps = 500, verbose = false)
    ok_La = maximum(sol_02.state.νt_c) > maximum(sol_03.state.νt_c)

    # KPPLT Langmuir 开关
    sol_kpp = run_to_steady(xuan_shen_config(; Nz = Nz, La_t = 0.2,
                                             closure = KPPLTClosure(; use_langmuir = true));
                            tol = 1e-10, verbose = false)
    sol_kpp0 = run_to_steady(xuan_shen_config(; Nz = Nz, La_t = 0.2,
                                              closure = KPPLTClosure(; use_langmuir = false));
                             tol = 1e-10, verbose = false)
    ok_kpp = maximum(sol_kpp.state.νt_c) > maximum(sol_kpp0.state.νt_c)

    ok = ok_E6 && ok_La && ok_kpp
    detail = @sprintf("νtmax E6=4/0: %.4e/%.4e; La0.2/0.3: %.4e/%.4e; KPPLT LT on/off: %.4e/%.4e",
                      ν_on, ν_off, maximum(sol_02.state.νt_c), maximum(sol_03.state.νt_c),
                      maximum(sol_kpp.state.νt_c), maximum(sol_kpp0.state.νt_c))
    return _pass("Langmuir enhances mixing (E6, La_t, KPPLT)", ok, detail;
                 metric = ν_on / max(ν_off, eps()))
end

# ---------------------------------------------------------------------------
# 6. KPPLT 形状函数峰值位置 σ=1/3
# ---------------------------------------------------------------------------
function check_kpp_shape_peak(; Nz = 96)
    cfg = xuan_shen_config(; Nz = Nz, La_t = 0.3,
                           closure = KPPLTClosure(; use_langmuir = true))
    sol = run_to_steady(cfg; tol = 1e-10, verbose = false)
    imax = argmax(sol.state.νt_c)
    σ = -cfg.grid.zc[imax] / cfg.grid.H
    ok = abs(σ - 1 / 3) < 2 / Nz + 0.02
    detail = @sprintf("νt peak at σ=%.4f (theory 1/3≈0.333)", σ)
    return _pass("KPPLT shape peak near σ=1/3", ok, detail; metric = abs(σ - 1 / 3))
end

# ---------------------------------------------------------------------------
# 7. 与 LES 数字化 νt 对比（形态趋势，非点对点强制相等）
# ---------------------------------------------------------------------------
function load_les_nut_csv(path::AbstractString)
    y = Float64[]
    n02 = Float64[]
    n03 = Float64[]
    open(path, "r") do io
        readline(io)  # header
        for line in eachline(io)
            isempty(strip(line)) && continue
            parts = split(strip(line), ',')
            push!(y, parse(Float64, parts[1]))
            push!(n02, parse(Float64, parts[2]))
            push!(n03, parse(Float64, parts[3]))
        end
    end
    return y, n02, n03
end

function _interp_linear(x::AbstractVector, f::AbstractVector, xq::Real)
    if xq >= x[1]
        return f[1]
    elseif xq <= x[end]
        return f[end]
    end
    # x is decreasing (0 → -1)
    for i in 1:length(x)-1
        if (x[i] >= xq >= x[i + 1]) || (x[i] <= xq <= x[i + 1])
            t = (xq - x[i]) / (x[i + 1] - x[i])
            return f[i] + t * (f[i + 1] - f[i])
        end
    end
    return f[end]
end

function pearson_corr(a::AbstractVector, b::AbstractVector)
    am = sum(a) / length(a)
    bm = sum(b) / length(b)
    num = sum((a[i] - am) * (b[i] - bm) for i in eachindex(a))
    den = sqrt(sum((a[i] - am)^2 for i in eachindex(a)) *
               sum((b[i] - bm)^2 for i in eachindex(b)))
    return den == 0 ? 0.0 : num / den
end

function check_les_nut_comparison(les_csv::AbstractString; Nz = 128, Cw = 0.55)
    y_les, n02_les, n03_les = load_les_nut_csv(les_csv)

    function rans_nut(La_t, closure)
        cfg = xuan_shen_config(; Nz = Nz, La_t = La_t, Reτ = 1000, k0H = 3.5, closure = closure)
        sol = run_to_steady(cfg; tol = 1e-7, verbose = false)
        y = cfg.grid.zc ./ cfg.grid.H
        ν = sol.state.νt_c ./ (cfg.forcing.u★ * cfg.grid.H)
        return y, ν, sol
    end

    results = CheckResult[]
    les_peak_02 = maximum(n02_les)
    les_peak_03 = maximum(n03_les)

    push!(results, _pass("LES reference loaded", true,
        @sprintf("N=%d  peak νt La0.2/0.3 = %.4f/%.4f", length(y_les), les_peak_02, les_peak_03);
        metric = les_peak_03))

    for (tag, clos) in (
            ("KLStokes", KLStokesClosure(; E6 = 4.0, channel = true)),
            ("KPPLT", KPPLTClosure(; Cw = Cw, use_langmuir = true)),
        )
        y02, ν02, _ = rans_nut(0.2, clos)
        y03, ν03, _ = rans_nut(0.3, clos)

        les02_i = [_interp_linear(y_les, n02_les, yi) for yi in y02]
        les03_i = [_interp_linear(y_les, n03_les, yi) for yi in y03]

        c02 = pearson_corr(ν02 ./ max(maximum(ν02), eps()),
                           les02_i ./ max(maximum(les02_i), eps()))
        c03 = pearson_corr(ν03 ./ max(maximum(ν03), eps()),
                           les03_i ./ max(maximum(les03_i), eps()))

        # 物理形态：单峰在内部；两端近壁涡粘显著小于峰值
        function structure_ok(y, ν)
            imax = argmax(ν)
            σ = -y[imax]
            peak = ν[imax]
            wall_ok = ν[1] < 0.25 * peak && ν[end] < 0.25 * peak
            interior_ok = 0.12 < σ < 0.88
            return wall_ok && interior_ok, σ, peak
        end
        ok02, σ02, p02 = structure_ok(y02, ν02)
        ok03, σ03, p03 = structure_ok(y03, ν03)

        # LES 对照：峰值深度误差（诊断）+ 结构通过为硬性
        iL = argmax(n03_les)
        dpeak = abs(y03[argmax(ν03)] - y_les[iL])
        ok = ok02 && ok03

        detail = @sprintf("%s: struct_ok=%s/%s  σ_peak(0.2/0.3)=%.3f/%.3f  corr=%.3f/%.3f  |Δσ_peak,LES|=%.3f  νtmax_R/L(0.3)=%.3f/%.3f",
                          tag, ok02, ok03, σ02, σ03, c02, c03, dpeak, p03, les_peak_03)
        push!(results, _pass("LES-informed νt structure ($tag)", ok, detail;
                             metric = min(c02, c03)))
    end

    return results
end

# ---------------------------------------------------------------------------
# 汇总运行
# ---------------------------------------------------------------------------
"""
    run_physics_validation(; les_csv, verbose=true) -> (passed, results)

运行全部物理验证项。任一硬性检查失败则 `passed=false`。
"""
function run_physics_validation(; les_csv::AbstractString = "",
                                verbose::Bool = true,
                                Nz::Int = 64)
    results = CheckResult[]

    cfg = xuan_shen_config(; Nz = Nz, La_t = 0.3, Reτ = 1000, E6 = 4.0)
    sol = run_to_steady(cfg; tol = 1e-7, max_steps = 500, verbose = false)
    push!(results, _pass("steady convergence (KLStokes La_t=0.3)", sol.converged,
                         @sprintf("iters=%d residual=%.3e", sol.iterations, sol.residual);
                         metric = sol.residual))

    push!(results, check_stress_balance(sol))
    push!(results, check_stokes_analytics(cfg))
    push!(results, check_velocity_structure(sol))
    push!(results, check_tke_production_dissipation(sol))
    push!(results, check_langmuir_mixing_trends(; Nz = min(Nz, 48)))
    push!(results, check_kpp_shape_peak(; Nz = max(Nz, 96)))

    if !isempty(les_csv) && isfile(les_csv)
        append!(results, check_les_nut_comparison(les_csv; Nz = max(Nz, 96)))
    else
        push!(results, _pass("LES νt comparison", true,
                             "skipped (LES CSV not provided)"; metric = 0.0))
    end

    if verbose
        println("="^72)
        println("Ocean1DRANS physics validation")
        println("="^72)
        for r in results
            mark = r.passed ? "PASS" : "FAIL"
            @printf("[%s] %s\n       %s\n", mark, r.name, r.detail)
        end
        npass = count(r -> r.passed, results)
        println("-"^72)
        @printf("%d / %d checks passed\n", npass, length(results))
    end

    passed = all(r -> r.passed, results)
    return passed, results
end
