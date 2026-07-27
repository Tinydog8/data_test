#!/usr/bin/env julia
#=
Xuan & Shen (2025) 稳态算例

论文 Fig.2(a) 是 Lagrangian 平均流 UL=U+Us。请对比 CSV 的 UL 与 Us。

闭合：
  :my25 / :kc04 — Mellor–Yamada 2.5 + Kantha–Clayson (2004)，默认成熟模型
  :harcourt     — 同上 + Lagrangian 动量应力
  :les          — Fig.2b 数字化 LES νt（论文对照）
  :kpplt        — 峰值校准 KPPLT
  :klstokes     — 简化代数 k–ℓ（仅趋势）
=#

using Ocean1DRANS
using Printf

outdir = joinpath(@__DIR__, "..", "output")
mkpath(outdir)

closures = (("my25", :my25), ("harcourt", :harcourt), ("les", :les),
            ("kpplt", :kpplt), ("klstokes", :klstokes))

for La_t in (0.2, 0.3)
    for (tag, clos) in closures
        @printf("\n========== La_t=%.1f  closure=%s ==========\n", La_t, tag)
        cfg = xuan_shen_config(; Nz = 96, La_t = La_t, Reτ = 1000, k0H = 3.5, closure = clos)
        sol = run_to_steady(cfg; tol = 1e-5, max_steps = 3000, verbose = true)
        csv = joinpath(outdir, @sprintf("xuan_shen_La%.1f_%s.csv", La_t, tag))
        write_profiles_csv(csv, sol)
        Umax = maximum(abs, sol.state.U)
        Usmax = maximum(abs, cfg.stokes.us_c)
        ULmax = maximum(abs, sol.state.U .+ cfg.stokes.us_c)
        @printf("Wrote %s\n", csv)
        @printf("max|U|=%.4f  max|Us|=%.4f  max|UL|=%.4f  max(νt)=%.4e  |U|/|Us|=%.3f\n",
                Umax, Usmax, ULmax, maximum(sol.state.νt_c), Umax / Usmax)
    end
end

@printf("\nDone. 成熟模型请看 my25/harcourt；论文对照看 les 的 UL 与 nu_t。\n")
