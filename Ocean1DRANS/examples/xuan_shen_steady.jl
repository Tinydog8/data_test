#!/usr/bin/env julia
#=
Xuan & Shen (2025) 稳态算例

论文 Fig.2(a) 是 Lagrangian 平均流 UL=U+Us。请对比 CSV 的 UL 与 Us。

闭合：
  :harcourt / :h15 — Harcourt (2015) 完整 SMC（默认，GOTM cmue_d_h15）
  :my25 / :kc04    — Mellor–Yamada 2.5 + Kantha–Clayson (2004)
  :les             — Fig.2b 数字化 LES νt（论文对照）
  :kpplt           — 峰值校准 KPPLT
  :klstokes        — 简化代数 k–ℓ（仅趋势）
=#

using Ocean1DRANS
using Printf

outdir = joinpath(@__DIR__, "..", "output")
mkpath(outdir)

closures = (("harcourt", :harcourt), ("my25", :my25), ("les", :les),
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
        imax = argmax(sol.state.νt_c)
        σ = -cfg.grid.zc[imax] / cfg.grid.H
        @printf("Wrote %s\n", csv)
        @printf("max|U|=%.4f  max|Us|=%.4f  max|UL|=%.4f  max(νt)=%.4e @σ=%.3f  max(νcl)=%.4e\n",
                Umax, Usmax, ULmax, maximum(sol.state.νt_c), σ, maximum(sol.state.νcl_c))
    end
end

@printf("\nDone. 推荐对照：harcourt（完整 SMC）与 les（Fig.2）；看 UL 与 nu_t。\n")
