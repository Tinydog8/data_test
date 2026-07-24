#!/usr/bin/env julia
#=
Xuan & Shen (2025) 型无分层 Langmuir 1D-RANS 稳态算例。

输出：
  - 稳态背景流 U(z), V(z)
  - 湍流粘性 ν_t(z)
  - Stokes 漂移 Us(z)
  - CSV 廓线文件，可供 resolvent / 后处理使用

运行（在 Ocean1DRANS 目录）：
  julia --project=. examples/xuan_shen_steady.jl
=#

using Ocean1DRANS
using Printf

outdir = joinpath(@__DIR__, "..", "output")
mkpath(outdir)

for La_t in (0.2, 0.3)
    @printf("\n========== La_t = %.1f (k–ℓ + Stokes production) ==========\n", La_t)
    cfg = xuan_shen_config(;
        Nz = 96,
        La_t = La_t,
        Reτ = 1000,
        k0H = 3.5,
        E6 = 4.0,
        nondimensional = true,
    )
    sol = run_to_steady(cfg; tol = 2e-5, max_steps = 300_000, check_every = 500, verbose = true)
    csv = joinpath(outdir, @sprintf("xuan_shen_La%.1f_klstokes.csv", La_t))
    write_profiles_csv(csv, sol)
    @printf("Wrote %s\n", csv)
    @printf("max|U|=%.4f  max(νt)=%.4e  Us(0)=%.4f\n",
            maximum(abs, sol.state.U), maximum(sol.state.νt_c), cfg.stokes.us_c[end])

    # 对照：KPPLT（一阶闭合）
    @printf("\n---------- La_t = %.1f (KPPLT) ----------\n", La_t)
    cfg_kpp = xuan_shen_config(;
        Nz = 96,
        La_t = La_t,
        Reτ = 1000,
        k0H = 3.5,
        nondimensional = true,
        closure = KPPLTClosure(; Cw = 0.15, use_langmuir = true),
    )
    sol_kpp = run_to_steady(cfg_kpp; tol = 2e-5, max_steps = 100_000,
                            check_every = 200, verbose = true)
    csv_kpp = joinpath(outdir, @sprintf("xuan_shen_La%.1f_kpplt.csv", La_t))
    write_profiles_csv(csv_kpp, sol_kpp)
    @printf("Wrote %s\n", csv_kpp)
end

@printf("\nDone. Profiles are in %s\n", outdir)
