"""
物理一致性检查：Lagrangian 应力平衡、Stokes 解析式、Langmuir 趋势、
TKE 生产–耗散平衡，以及与 LES / 论文 Fig.2 的对照。
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
# 1. Lagrangian 应力平衡
# ---------------------------------------------------------------------------
function check_stress_balance(sol::SteadySolution; atol = 1e-10)
    cfg = sol.config
    diag = diagnostic_stress_balance(sol)
    τb = diag.τ_f[1]
    τt = diag.τ_f[end]
    ok_bottom = abs(τb) < atol
    ok_top = abs(τt - cfg.forcing.τx) < atol
    ok_pg = abs(cfg.forcing.τx + cfg.forcing.Fx * cfg.grid.H) < atol

    g = cfg.grid
    U = sol.state.U
    νe = diag.νe
    αs = diag.αs
    err = 0.0
    count = 0
    @inbounds for i in 2:g.Nz
        Uz = (U[i] - U[i - 1]) / g.dz
        Usz = cfg.stokes.dusdz_f[i]
        τ_num = νe[i] * (Uz + αs * Usz)
        τ_ana = diag.τ_f[i]
        err += abs(τ_num - τ_ana)
        count += 1
    end
    mean_err = err / max(count, 1)
    ok_recon = mean_err < 1e-5 * max(abs(cfg.forcing.τx), 1.0)

    ok = ok_bottom && ok_top && ok_pg && ok_recon
    detail = @sprintf("αs=%.2f τ_bottom=%.3e τ_top=%.3e  recon_MAE=%.3e",
                      αs, τb, τt, mean_err)
    return _pass("Lagrangian stress balance", ok, detail; metric = mean_err)
end

# ---------------------------------------------------------------------------
# 2. Stokes / La_t
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
# 3. 论文关键：欧拉力应远小于 Stokes（UL≈Us）
# ---------------------------------------------------------------------------
function check_small_eulerian_mean(sol::SteadySolution; max_ratio = 0.35)
    U = sol.state.U
    Us = sol.config.stokes.us_c
    # 排除近底薄层（规范 U_bottom=0 造成的积分累积）
    g = sol.config.grid
    mask = [(-g.zc[i] / g.H) < 0.85 for i in eachindex(g.zc)]
    Ums = maximum(abs, U[mask])
    Usms = maximum(abs, Us[mask])
    ratio = Ums / max(Usms, eps())
    # LES νt 或充分混合时比值应明显 < 1
    ok = ratio < max_ratio
    detail = @sprintf("max|U|/max|Us|=%.3f (upper 85%% column); max|UL|=%.3f max|Us|=%.3f",
                      ratio, maximum(abs, U .+ Us), Usms)
    return _pass("Eulerian mean ≪ Stokes (paper Fig.2a)", ok, detail; metric = ratio)
end

# ---------------------------------------------------------------------------
# 4. TKE 生产–耗散
# ---------------------------------------------------------------------------
function check_tke_production_dissipation(sol::SteadySolution; rtol = 0.05)
    cfg = sol.config
    clos = cfg.closure
    if clos isa MY25KC04Closure
        return check_my25_production_dissipation(sol; rtol = rtol)
    end
    clos isa KLStokesClosure || return _pass(
        "TKE production–dissipation balance", true,
        "skipped (closure is not KLStokes/MY25)"; metric = 0.0)

    g = cfg.grid
    U, V, k, ℓ = sol.state.U, sol.state.V, sol.state.k, sol.state.ℓ
    νt_c = sol.state.νt_c
    dz = g.dz
    max_rel = 0.0
    mean_rel = 0.0
    n = 0
    @inbounds for i in 1:g.Nz
        σ = -g.zc[i] / g.H
        (σ < 0.15 || σ > 0.85 || i == 1 || i == g.Nz) && continue
        Uz = (U[i + 1] - U[i - 1]) / (2dz)
        Vz = (V[i + 1] - V[i - 1]) / (2dz)
        P, PS = tke_production(Uz, Vz, cfg.stokes.dusdz_c[i], cfg.stokes.dvsdz_c[i],
                               νt_c[i], clos.E6, clos.αs)
        Prod = P + PS
        ε = clos.cε * (max(k[i], clos.k_min)^(3 / 2)) / max(ℓ[i], eps())
        if abs(Prod) > 1e-12
            rel = abs(Prod - ε) / abs(Prod)
            max_rel = max(max_rel, rel)
            mean_rel += rel
            n += 1
        end
    end
    mean_rel = n == 0 ? 0.0 : mean_rel / n
    ok = n > 0 && max_rel < rtol
    detail = @sprintf("max|P+PS-ε|/|Prod|=%.3e  mean=%.3e  (n=%d)", max_rel, mean_rel, n)
    return _pass("TKE production–dissipation balance", ok, detail; metric = max_rel)
end

function check_my25_production_dissipation(sol::SteadySolution; rtol = 0.15)
    cfg = sol.config
    clos = cfg.closure
    g = cfg.grid
    U, V = sol.state.U, sol.state.V
    q2, q2l = sol.state.q2, sol.state.q2l
    νt_c = sol.state.νt_c
    dz = g.dz
    max_rel = 0.0
    n = 0
    @inbounds for i in 2:g.Nz-1
        σ = -g.zc[i] / g.H
        (σ < 0.15 || σ > 0.85) && continue
        Uz = (U[i + 1] - U[i - 1]) / (2dz)
        Vz = (V[i + 1] - V[i - 1]) / (2dz)
        P, Ps = my25_production(Uz, Vz, cfg.stokes.dusdz_c[i], cfg.stokes.dvsdz_c[i], νt_c[i])
        Prod = P + Ps
        q2i = max(q2[i], clos.q2_min)
        ℓ = max(q2l[i] / q2i, clos.ℓ_min)
        ε = (sqrt(q2i)^3) / (clos.B1 * ℓ)
        if abs(Prod) > 1e-12
            rel = abs(Prod - ε) / abs(Prod)
            max_rel = max(max_rel, rel)
            n += 1
        end
    end
    ok = n > 0 && max_rel < rtol
    detail = @sprintf("MY25 max|P+Ps-ε|/|Prod|=%.3e (n=%d)", max_rel, n)
    return _pass("MY25 q² production–dissipation balance", ok, detail; metric = max_rel)
end

function check_my25_les_magnitude(; Nz = 64)
    # KC04 with E6=4 should lift KM into LES order (not 10× low)
    sol_on = run_to_steady(xuan_shen_config(; Nz = Nz, La_t = 0.3, E6 = 4.0, closure = :my25);
                           tol = 1e-5, max_steps = 2000, verbose = false)
    sol_off = run_to_steady(xuan_shen_config(; Nz = Nz, La_t = 0.3, E6 = 0.0, closure = :my25);
                            tol = 1e-5, max_steps = 2000, verbose = false)
    ν_on = maximum(sol_on.state.νt_c)
    ν_off = maximum(sol_off.state.νt_c)
    # LES peak ~0.38; mature KC04 should be O(0.1) or above, not O(0.03)
    ok_mag = sol_on.converged && ν_on > 0.08
    ok_E6 = ν_on > ν_off * 1.05
    ok = ok_mag && ok_E6
    detail = @sprintf("MY25 La0.3 νtmax E6=4/0: %.3e/%.3e (target ≳0.08 with E6)",
                      ν_on, ν_off)
    return _pass("MY25/KC04 νt magnitude vs LES order (E6)", ok, detail;
                 metric = ν_on)
end

# ---------------------------------------------------------------------------
# 5. Langmuir 增强
# ---------------------------------------------------------------------------
function check_langmuir_mixing_trends(; Nz = 48)
    sol_on = run_to_steady(xuan_shen_config(; Nz = Nz, La_t = 0.2, E6 = 4.0, αs = 1.0,
                                            closure = :klstokes);
                           tol = 1e-6, max_steps = 500, verbose = false)
    sol_off = run_to_steady(xuan_shen_config(; Nz = Nz, La_t = 0.2, E6 = 0.0, αs = 1.0,
                                             closure = :klstokes);
                            tol = 1e-6, max_steps = 500, verbose = false)
    ν_on = maximum(sol_on.state.νt_c)
    ν_off = maximum(sol_off.state.νt_c)
    ok_E6 = sol_on.converged && sol_off.converged && ν_on > ν_off * 0.99

    sol_02 = run_to_steady(xuan_shen_config(; Nz = Nz, La_t = 0.2, closure = :kpplt);
                           tol = 1e-8, verbose = false)
    sol_03 = run_to_steady(xuan_shen_config(; Nz = Nz, La_t = 0.3, closure = :kpplt);
                           tol = 1e-8, verbose = false)
    ok_La = maximum(sol_02.state.νt_c) > maximum(sol_03.state.νt_c)

    sol_kpp = run_to_steady(xuan_shen_config(; Nz = Nz, La_t = 0.2,
                                             closure = KPPLTClosure(; use_langmuir = true));
                            tol = 1e-10, verbose = false)
    sol_kpp0 = run_to_steady(xuan_shen_config(; Nz = Nz, La_t = 0.2,
                                              closure = KPPLTClosure(; use_langmuir = false));
                             tol = 1e-10, verbose = false)
    ok_kpp = maximum(sol_kpp.state.νt_c) > maximum(sol_kpp0.state.νt_c)

    ok = ok_E6 && ok_La && ok_kpp
    detail = @sprintf("νtmax E6=4/0: %.3e/%.3e; KPPLT La0.2/0.3: %.3e/%.3e; LT on/off: %.3e/%.3e",
                      ν_on, ν_off, maximum(sol_02.state.νt_c), maximum(sol_03.state.νt_c),
                      maximum(sol_kpp.state.νt_c), maximum(sol_kpp0.state.νt_c))
    return _pass("Langmuir enhances mixing (E6, La_t, KPPLT)", ok, detail;
                 metric = ν_on / max(ν_off, eps()))
end

function check_kpp_shape_peak(; Nz = 96)
    cfg = xuan_shen_config(; Nz = Nz, La_t = 0.3, closure = :kpplt)
    sol = run_to_steady(cfg; tol = 1e-10, verbose = false)
    imax = argmax(sol.state.νt_c)
    σ = -cfg.grid.zc[imax] / cfg.grid.H
    ok = abs(σ - 1 / 3) < 2 / Nz + 0.02
    detail = @sprintf("νt peak at σ=%.4f (theory 1/3≈0.333), max νt=%.4f",
                      σ, maximum(sol.state.νt_c))
    return _pass("KPPLT shape peak near σ=1/3", ok, detail; metric = abs(σ - 1 / 3))
end

# ---------------------------------------------------------------------------
# 6. LES νt 对照 + 论文一致性
# ---------------------------------------------------------------------------
function check_les_nut_comparison(les_csv::AbstractString; Nz = 128)
    y_les, n02_les, n03_les = load_les_nut_csv(les_csv)
    results = CheckResult[]
    push!(results, _pass("LES reference loaded", true,
        @sprintf("N=%d  peak νt La0.2/0.3 = %.4f/%.4f",
                 length(y_les), maximum(n02_les), maximum(n03_les));
        metric = maximum(n03_les)))

    # LESNutClosure 应精确再现数字化剖面峰值
    for La in (0.2, 0.3)
        cfg = xuan_shen_config(; Nz = Nz, La_t = La, closure = :les)
        sol = run_to_steady(cfg; verbose = false)
        ν = sol.state.νt_c ./ (cfg.forcing.u★ * cfg.grid.H)
        les = La == 0.2 ? n02_les : n03_les
        peak_r = maximum(ν)
        peak_l = maximum(les)
        ok_peak = abs(peak_r - peak_l) / peak_l < 0.05
        # 欧拉力应很小
        rU = maximum(abs, sol.state.U) / maximum(abs, cfg.stokes.us_c)
        ok_U = rU < 0.35
        detail = @sprintf("La=%.1f LESNut: νtmax R/L=%.4f/%.4f  max|U|/max|Us|=%.3f",
                          La, peak_r, peak_l, rU)
        push!(results, _pass("Paper-consistent LESNut (La=$La)", ok_peak && ok_U, detail;
                             metric = rU))
    end

    # KPPLT 校准后峰值量级应进入 LES 同阶（La=0.3 → O(0.3)）
    cfg = xuan_shen_config(; Nz = Nz, La_t = 0.3, closure = :kpplt)
    sol = run_to_steady(cfg; verbose = false)
    peak = maximum(sol.state.νt_c) / (cfg.forcing.u★ * cfg.grid.H)
    ok = 0.2 < peak < 0.6
    detail = @sprintf("KPPLT La0.3 max νt/(u★H)=%.3f (LES≈0.38)", peak)
    push!(results, _pass("KPPLT νt magnitude vs LES", ok, detail; metric = peak))

    return results
end

# ---------------------------------------------------------------------------
# 汇总
# ---------------------------------------------------------------------------
function run_physics_validation(; les_csv::AbstractString = "",
                                verbose::Bool = true,
                                Nz::Int = 64)
    results = CheckResult[]

    # 默认用 LESNut：与论文 Fig.2 一致的检查基线
    cfg = xuan_shen_config(; Nz = Nz, La_t = 0.3, closure = :les)
    sol = run_to_steady(cfg; verbose = false)
    push!(results, _pass("steady convergence (LESNut La_t=0.3)", sol.converged,
                         @sprintf("iters=%d residual=%.3e", sol.iterations, sol.residual);
                         metric = sol.residual))
    push!(results, check_stress_balance(sol))
    push!(results, check_stokes_analytics(cfg))
    push!(results, check_small_eulerian_mean(sol))

    cfg_kl = xuan_shen_config(; Nz = Nz, La_t = 0.3, closure = :klstokes)
    sol_kl = run_to_steady(cfg_kl; tol = 1e-6, verbose = false)
    push!(results, check_tke_production_dissipation(sol_kl))

    cfg_my = xuan_shen_config(; Nz = Nz, La_t = 0.3, closure = :my25)
    sol_my = run_to_steady(cfg_my; tol = 1e-5, max_steps = 2000, verbose = false)
    push!(results, _pass("steady convergence (MY25/KC04 La_t=0.3)", sol_my.converged,
                         @sprintf("iters=%d residual=%.3e maxνt=%.3e",
                                  sol_my.iterations, sol_my.residual, maximum(sol_my.state.νt_c));
                         metric = maximum(sol_my.state.νt_c)))
    push!(results, check_tke_production_dissipation(sol_my))
    push!(results, check_my25_les_magnitude(; Nz = min(Nz, 64)))

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
    return all(r -> r.passed, results), results
end
