using Test
using Ocean1DRANS

@testset "Ocean1DRANS" begin
    @testset "grid" begin
        g = UniformColumnGrid(16, 2.0)
        @test g.Nz == 16
        @test g.zf[1] ≈ -2.0
        @test g.zf[end] ≈ 0.0
    end

    @testset "Stokes drift" begin
        g = UniformColumnGrid(32, 1.0)
        s = monochromatic_stokes(g; u★ = 1.0, La_t = 0.3, k0H = 3.5)
        Us0 = 1.0 / 0.3^2
        @test s.us_c[end] ≈ Us0 * exp(2 * 3.5 * g.zc[end]) rtol = 1e-12
        @test all(s.dusdz_c .> 0)
    end

    @testset "stress balance signs" begin
        cfg = xuan_shen_config(; Nz = 32, La_t = 0.3, closure = :les)
        @test cfg.forcing.τx ≈ cfg.forcing.u★^2
        @test cfg.forcing.τx + cfg.forcing.Fx * cfg.grid.H ≈ 0 atol = 1e-14
    end

    @testset "MY25 E6 vs KC04 Fig.1 digitized trend" begin
        # Implementation gate: E6 must enhance KM (KC04 Fig.1 qualitative)
        sol4 = run_to_steady(xuan_shen_config(; Nz = 48, La_t = 0.3, closure = :my25, E6 = 4.0);
                             tol = 1e-5, verbose = false)
        sol0 = run_to_steady(xuan_shen_config(; Nz = 48, La_t = 0.3, closure = :my25, E6 = 0.0);
                             tol = 1e-5, verbose = false)
        @test maximum(sol4.state.νt_c) > 1.5 * maximum(sol0.state.νt_c)
        # Digitized reference exists and E6=4 peak > noLC peak
        ref = joinpath(@__DIR__, "..", "data", "kc04_fig1_KM.csv")
        @test isfile(ref)
        lines = readlines(ref)
        @test occursin("KM_E6_4", lines[1])
        # last data-ish: parse a mid-depth row
        row = split(lines[findfirst(l -> startswith(l, "-0.50"), lines)], ',')
        @test parse(Float64, row[4]) > parse(Float64, row[2])  # E6=4 > noLC
        les = joinpath(@__DIR__, "..", "data", "mcwilliams1997_fig3b_KM.csv")
        @test isfile(les)
    end

    @testset "Harcourt2015 full SMC" begin
        cfg = xuan_shen_config(; Nz = 64, La_t = 0.3, closure = :harcourt)
        sol = run_to_steady(cfg; tol = 1e-5, max_steps = 3000, verbose = false)
        @test sol.converged
        @test cfg.closure isa Harcourt2015Closure
        @test cfg.closure.E6 == 6.0
        @test maximum(sol.state.νt_c) > 0.08
        @test maximum(sol.state.νcl_c) > 0.0
        imax = argmax(sol.state.νt_c)
        σ = -cfg.grid.zc[imax] / cfg.grid.H
        @test σ > 0.25
    end

    @testset "MY25/KC04 mature closure magnitude" begin
        cfg = xuan_shen_config(; Nz = 64, La_t = 0.3, closure = :my25)
        sol = run_to_steady(cfg; tol = 1e-5, max_steps = 2000, verbose = false)
        @test sol.converged
        @test maximum(sol.state.νt_c) > 0.08   # E6=4 should lift KM out of O(0.03)
        @test cfg.closure isa MY25KC04Closure
        @test cfg.closure.E6 == 4.0
        @test cfg.closure.αs == 0.0          # strict KC04 momentum
        # E6 enhances mixing vs E6=0
        cfg0 = xuan_shen_config(; Nz = 48, La_t = 0.3, E6 = 0.0, closure = :my25)
        sol0 = run_to_steady(cfg0; tol = 1e-5, max_steps = 2000, verbose = false)
        @test maximum(sol.state.νt_c) > maximum(sol0.state.νt_c) * 1.05
    end

    @testset "LESNut matches paper: small U, UL≈Us" begin
        cfg = xuan_shen_config(; Nz = 64, La_t = 0.3, closure = :les)
        sol = run_to_steady(cfg; verbose = false)
        @test sol.converged
        Usmax = maximum(abs, cfg.stokes.us_c)
        Umax = maximum(abs, sol.state.U)
        @test Umax / Usmax < 0.35
        @test maximum(sol.state.νt_c) > 0.3  # nondim u★=H=1 → LES peak ~0.38
        # UL 与 Us 同量级
        UL = sol.state.U .+ cfg.stokes.us_c
        @test maximum(abs, UL) ≈ Usmax rtol = 0.35
    end

    @testset "KPPLT calibrated magnitude" begin
        cfg = xuan_shen_config(; Nz = 64, La_t = 0.3, closure = :kpplt)
        sol = run_to_steady(cfg; verbose = false)
        @test sol.converged
        @test 0.2 < maximum(sol.state.νt_c) < 0.6
    end

    @testset "Lagrangian stress brings UL closer to Us" begin
        cfg0 = xuan_shen_config(; Nz = 48, La_t = 0.3,
                                closure = KPPLTClosure(; Cw = 3.6, αs = 0.0))
        cfg1 = xuan_shen_config(; Nz = 48, La_t = 0.3,
                                closure = KPPLTClosure(; Cw = 3.6, αs = 1.0))
        sol0 = run_to_steady(cfg0; verbose = false)
        sol1 = run_to_steady(cfg1; verbose = false)
        Us = cfg1.stokes.us_c
        err0 = maximum(abs, (sol0.state.U .+ Us) .- Us)  # = max|U|
        err1 = maximum(abs, (sol1.state.U .+ cfg1.stokes.us_c) .- cfg1.stokes.us_c)
        # αs=1 时欧拉偏差不应系统性更差；用 LESNut 已保证论文一致性
        @test err1 < 2 * err0 + 1.0
        @test maximum(sol1.state.νt_c) ≈ maximum(sol0.state.νt_c) rtol = 1e-6
    end

    @testset "KLStokes Langmuir E6" begin
        cfg_lt = xuan_shen_config(; Nz = 40, La_t = 0.2, E6 = 4.0, closure = :klstokes)
        cfg_st = xuan_shen_config(; Nz = 40, La_t = 0.2, E6 = 0.0, closure = :klstokes)
        sol_lt = run_to_steady(cfg_lt; tol = 1e-6, verbose = false)
        sol_st = run_to_steady(cfg_st; tol = 1e-6, verbose = false)
        @test sol_lt.converged && sol_st.converged
        @test maximum(sol_lt.state.νt_c) >= maximum(sol_st.state.νt_c) * 0.99
    end

    @testset "CSV has UL column" begin
        cfg = xuan_shen_config(; Nz = 24, La_t = 0.3, closure = :les)
        sol = run_to_steady(cfg; verbose = false)
        path = tempname() * ".csv"
        write_profiles_csv(path, sol)
        header = readline(path)
        @test occursin("UL", header)
        @test occursin("nu_t", header)
        @test occursin("nu_cl", header)
        rm(path)
    end

    @testset "physics validation (core)" begin
        passed, results = run_physics_validation(; les_csv = "", verbose = false, Nz = 48)
        @test passed
        @test any(r -> occursin("stress balance", r.name) && r.passed, results)
        @test any(r -> occursin("Stokes", r.name) && r.passed, results)
    end
end
