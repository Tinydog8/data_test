#!/usr/bin/env julia
#=
物理准确性验证（含与 Xuan & Shen 2025 Fig.2 的对照）

运行：
  julia --project=. examples/validate_physics.jl
=#

using Ocean1DRANS
using Printf

les_csv = joinpath(@__DIR__, "..", "data", "les_eddy_viscosity_fig2b.csv")
outdir = joinpath(@__DIR__, "..", "output")
mkpath(outdir)

println("LES reference: ", les_csv, "  exists=", isfile(les_csv))
passed, results = run_physics_validation(; les_csv = les_csv, verbose = true, Nz = 64)

# 写出论文对照廓线：UL vs Us，νt vs LES / Harcourt / MY25 / KPPLT
function write_paper_comparison(path, La_t)
    y_les, n02, n03 = Ocean1DRANS.load_les_nut_csv(les_csv)
    n_les = La_t ≈ 0.2 ? n02 : n03
    open(path, "w") do io
        println(io, "y_over_H,UL_les_proxy_Us,UL_RANS_lesNut,UL_RANS_harcourt,UL_RANS_my25,UL_RANS_kpplt,nu_t_LES,nu_t_lesNut,nu_t_harcourt,nu_t_my25,nu_t_kpplt,U_euler_harcourt")
        cfg_les = xuan_shen_config(; Nz = 128, La_t = La_t, closure = :les)
        sol_les = run_to_steady(cfg_les; verbose = false)
        cfg_h = xuan_shen_config(; Nz = 128, La_t = La_t, closure = :harcourt)
        sol_h = run_to_steady(cfg_h; tol = 1e-5, max_steps = 3000, verbose = false)
        cfg_my = xuan_shen_config(; Nz = 128, La_t = La_t, closure = :my25)
        sol_my = run_to_steady(cfg_my; tol = 1e-5, max_steps = 3000, verbose = false)
        cfg_kpp = xuan_shen_config(; Nz = 128, La_t = La_t, closure = :kpplt)
        sol_kpp = run_to_steady(cfg_kpp; verbose = false)
        uH = cfg_les.forcing.u★ * cfg_les.grid.H
        for i in eachindex(cfg_les.grid.zc)
            yH = cfg_les.grid.zc[i] / cfg_les.grid.H
            Us = cfg_les.stokes.us_c[i]
            UL_les = sol_les.state.U[i] + Us
            UL_h = sol_h.state.U[i] + cfg_h.stokes.us_c[i]
            UL_my = sol_my.state.U[i] + cfg_my.stokes.us_c[i]
            UL_kpp = sol_kpp.state.U[i] + cfg_kpp.stokes.us_c[i]
            ν_lesn = sol_les.state.νt_c[i] / uH
            ν_h = sol_h.state.νt_c[i] / (cfg_h.forcing.u★ * cfg_h.grid.H)
            ν_my = sol_my.state.νt_c[i] / (cfg_my.forcing.u★ * cfg_my.grid.H)
            ν_kpp = sol_kpp.state.νt_c[i] / (cfg_kpp.forcing.u★ * cfg_kpp.grid.H)
            ν_ref = Ocean1DRANS._interp_linear(y_les, n_les, yH)
            @printf(io, "%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e\n",
                    yH, Us, UL_les, UL_h, UL_my, UL_kpp, ν_ref, ν_lesn, ν_h, ν_my, ν_kpp, sol_h.state.U[i])
        end
    end
    @printf("Wrote %s\n", path)
end

if isfile(les_csv)
    write_paper_comparison(joinpath(outdir, "paper_compare_La0.2.csv"), 0.2)
    write_paper_comparison(joinpath(outdir, "paper_compare_La0.3.csv"), 0.3)
end

summary_path = joinpath(outdir, "validation_summary.txt")
open(summary_path, "w") do io
    println(io, "Ocean1DRANS physics validation summary")
    println(io, "passed=$(passed)")
    for r in results
        @printf(io, "%s\t%s\t%.6e\t%s\n", r.passed ? "PASS" : "FAIL", r.name, r.metric, r.detail)
    end
end
@printf("Wrote %s\n", summary_path)

if !passed
    println("\nVALIDATION FAILED")
    exit(1)
end
println("\nVALIDATION PASSED")
