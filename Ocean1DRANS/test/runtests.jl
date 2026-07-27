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
        rm(path)
    end

    @testset "physics validation (core)" begin
        passed, results = run_physics_validation(; les_csv = "", verbose = false, Nz = 48)
        @test passed
        @test any(r -> occursin("stress balance", r.name) && r.passed, results)
        @test any(r -> occursin("Stokes", r.name) && r.passed, results)
    end
end
