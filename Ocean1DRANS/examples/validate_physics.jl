#!/usr/bin/env julia
#=
物理准确性验证脚本
==================

检查 Ocean1DRANS 是否满足关键物理约束，并与 Xuan & Shen (2025)
Fig.2(b) LES 数字化涡粘剖面做形态对比。

硬性检查（失败则 exit 1）：
  1. 稳态收敛
  2. 动量应力平衡（底应力 0、表应力 = u★²、离散重建）
  3. Stokes 指数廓线与 La_t = √(u★/Us0) 定义
  4. 风生剪切结构（∂U/∂z ≥ 0）
  5. k–ℓ 局部生产–耗散平衡 P + E6·P_S ≈ ε
  6. Langmuir 增强混合（E6、La_t、KPPLT 开关）
  7. KPPLT 形状峰位于 σ ≈ 1/3
  8. 与 LES νt 归一化形状相关（Pearson > 0.7）

运行：
  cd Ocean1DRANS
  julia --project=. examples/validate_physics.jl
=#

using Ocean1DRANS
using Printf

les_csv = joinpath(@__DIR__, "..", "data", "les_eddy_viscosity_fig2b.csv")
outdir = joinpath(@__DIR__, "..", "output")
mkpath(outdir)

println("LES reference: ", les_csv, "  exists=", isfile(les_csv))

passed, results = run_physics_validation(; les_csv = les_csv, verbose = true, Nz = 64)

# 额外写出对照廓线，便于目视检查
function write_comparison_csv(path, La_t)
    y_les, n02, n03 = Ocean1DRANS.load_les_nut_csv(les_csv)
    n_les = La_t ≈ 0.2 ? n02 : n03

    cfg_kl = xuan_shen_config(; Nz = 128, La_t = La_t, E6 = 4.0)
    sol_kl = run_to_steady(cfg_kl; tol = 1e-7, verbose = false)
    cfg_kpp = xuan_shen_config(; Nz = 128, La_t = La_t,
                               closure = KPPLTClosure(; Cw = 0.55))
    sol_kpp = run_to_steady(cfg_kpp; tol = 1e-10, verbose = false)

    open(path, "w") do io
        println(io, "y_over_H,nu_t_LES,nu_t_KLStokes,nu_t_KPPLT")
        for i in eachindex(cfg_kl.grid.zc)
            yH = cfg_kl.grid.zc[i] / cfg_kl.grid.H
            ν_kl = sol_kl.state.νt_c[i] / (cfg_kl.forcing.u★ * cfg_kl.grid.H)
            ν_kpp = sol_kpp.state.νt_c[i] / (cfg_kpp.forcing.u★ * cfg_kpp.grid.H)
            ν_les = Ocean1DRANS._interp_linear(y_les, n_les, yH)
            @printf(io, "%.8e,%.8e,%.8e,%.8e\n", yH, ν_les, ν_kl, ν_kpp)
        end
    end
    @printf("Wrote comparison profiles: %s\n", path)
end

if isfile(les_csv)
    write_comparison_csv(joinpath(outdir, "validate_nut_La0.2.csv"), 0.2)
    write_comparison_csv(joinpath(outdir, "validate_nut_La0.3.csv"), 0.3)
end

# 机器可读摘要
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
