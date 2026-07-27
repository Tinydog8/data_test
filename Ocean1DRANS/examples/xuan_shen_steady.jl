#!/usr/bin/env julia
#=
Xuan & Shen (2025) 稳态算例（修正版）

重要：论文 Fig.2(a) 画的是 Lagrangian 平均流 UL=U+Us，且欧拉力几乎可忽略。
请对比 CSV 中的 **UL** 列与 Us，而不是欧拉 U。

三种闭合：
  :les      — Fig.2b 数字化 LES νt（与论文最一致，推荐作 resolvent 基流）
  :kpplt    — 峰值校准的 KPPLT（αs=1 Lagrangian 应力）
  :klstokes — k–ℓ + Stokes 生产 + Lagrangian 应力
=#

using Ocean1DRANS
using Printf

outdir = joinpath(@__DIR__, "..", "output")
mkpath(outdir)

for La_t in (0.2, 0.3)
    for (tag, clos) in (("les", :les), ("kpplt", :kpplt), ("klstokes", :klstokes))
        @printf("\n========== La_t=%.1f  closure=%s ==========\n", La_t, tag)
        cfg = xuan_shen_config(; Nz = 96, La_t = La_t, Reτ = 1000, k0H = 3.5, closure = clos)
        sol = run_to_steady(cfg; tol = 1e-6, verbose = true)
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

@printf("\nDone. 对照论文时请使用 UL 与 nu_t 列。\n")
