using Test
using Ocean1DRANS

@testset "Ocean1DRANS" begin
    @testset "grid" begin
        g = UniformColumnGrid(16, 2.0)
        @test g.Nz == 16
        @test g.H == 2.0
        @test length(g.zc) == 16
        @test length(g.zf) == 17
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
        cfg = xuan_shen_config(; Nz = 32, La_t = 0.3)
        @test cfg.forcing.τx ≈ cfg.forcing.u★^2
        @test cfg.forcing.Fx ≈ -cfg.forcing.u★^2 / cfg.grid.H
        # 底应力应为 0：τx - Fx*(-H) = τx + Fx*H = 0
        @test cfg.forcing.τx + cfg.forcing.Fx * cfg.grid.H ≈ 0 atol = 1e-14
    end

    @testset "Xuan-Shen steady KLStokes" begin
        cfg = xuan_shen_config(; Nz = 48, La_t = 0.3, Reτ = 1000, E6 = 4.0)
        sol = run_to_steady(cfg; tol = 1e-6, max_steps = 500, verbose = false)
        @test sol.converged
        @test maximum(sol.state.νt_c) > 0
        @test all(isfinite, sol.state.U)
        @test all(sol.state.k .> 0)
        @test cfg.stokes.us_f[end] ≈ 1 / 0.3^2 rtol = 1e-12
        imax = argmax(sol.state.νt_c)
        @test 4 < imax < cfg.grid.Nz - 3
        diag = Ocean1DRANS.diagnostic_stress_balance(sol)
        @test diag.τ_f[1] ≈ 0 atol = 1e-12
        @test diag.τ_f[end] ≈ 1 atol = 1e-12
    end

    @testset "KPPLT steady" begin
        cfg = xuan_shen_config(;
            Nz = 48, La_t = 0.2, Reτ = 1000,
            closure = KPPLTClosure(; use_langmuir = true),
        )
        sol = run_to_steady(cfg; tol = 1e-10, max_steps = 50, verbose = false)
        @test sol.converged
        @test maximum(sol.state.νt_c) > 0

        cfg0 = xuan_shen_config(;
            Nz = 48, La_t = 0.2, Reτ = 1000,
            closure = KPPLTClosure(; use_langmuir = false),
        )
        sol0 = run_to_steady(cfg0; tol = 1e-10, max_steps = 50, verbose = false)
        @test maximum(sol.state.νt_c) > maximum(sol0.state.νt_c)
    end

    @testset "Langmuir increases mixing (KLStokes E6)" begin
        cfg_lt = xuan_shen_config(; Nz = 40, La_t = 0.2, E6 = 4.0)
        cfg_st = xuan_shen_config(; Nz = 40, La_t = 0.2, E6 = 0.0)
        sol_lt = run_to_steady(cfg_lt; tol = 1e-6, max_steps = 500, verbose = false)
        sol_st = run_to_steady(cfg_st; tol = 1e-6, max_steps = 500, verbose = false)
        @test sol_lt.converged && sol_st.converged
        @test maximum(sol_lt.state.νt_c) > maximum(sol_st.state.νt_c)
    end

    @testset "CSV IO" begin
        cfg = xuan_shen_config(; Nz = 24, La_t = 0.3)
        sol = run_to_steady(cfg; tol = 1e-6, max_steps = 500, verbose = false)
        path = tempname() * ".csv"
        write_profiles_csv(path, sol)
        @test isfile(path)
        lines = readlines(path)
        @test startswith(lines[1], "z,U,V")
        @test length(lines) == 25
        rm(path)
    end
end
