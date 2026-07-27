#!/usr/bin/env julia
#=
对照 Kantha & Clayson (2004) Fig.1 / McWilliams et al. (1997) Fig.3b
验证 MY2.5 + KC04 的 E6 Langmuir 项。

硬指标（实现是否正确）：
  Xuan–Shen 通道应力平衡上，E6=4 必须明显高于 E6=0（ratio > 1.5）。

不要用这些判断 E6：
  未充分 spin-up / 无分层的 McWilliams+Coriolis 跑次——那里 Uz≈0 ⇒ Ps≈0，
  E6 几乎无效（旧脚本叠画列 KM_MY25_E6_* 就踩过这个坑）。

运行：
  julia --project=. examples/compare_my25_kc04.jl
=#

using Ocean1DRANS
using Printf

root = joinpath(@__DIR__, "..")
outdir = joinpath(root, "output")
mkpath(outdir)

function load_km_csv(path; zcol = 1, cols = 2:4)
    zs = Float64[]
    cols_data = [Float64[] for _ in cols]
    open(path) do io
        readline(io)
        for line in eachline(io)
            isempty(strip(line)) && continue
            p = split(strip(line), ',')
            push!(zs, parse(Float64, p[zcol]))
            for (j, c) in enumerate(cols)
                push!(cols_data[j], parse(Float64, p[c]))
            end
        end
    end
    return zs, cols_data
end

function interp1(x, f, xq)
    if xq >= maximum(x)
        return x[1] >= x[end] ? f[1] : f[end]
    elseif xq <= minimum(x)
        return x[1] >= x[end] ? f[end] : f[1]
    end
    for i in 1:length(x)-1
        x1, x2 = x[i], x[i + 1]
        if (x1 >= xq >= x2) || (x1 <= xq <= x2)
            t = (xq - x1) / (x2 - x1)
            return f[i] + t * (f[i + 1] - f[i])
        end
    end
    return f[end]
end

# ---------------------------------------------------------------------------
# A) 通道：E6 硬指标（叠画 MY25 列也来自这里）
# ---------------------------------------------------------------------------
println("="^72)
println("A) Channel MY25 E6=4 vs E6=0  [implementation gate]")
println("="^72)

sol4 = run_to_steady(xuan_shen_config(; Nz = 96, La_t = 0.3, closure = :my25, E6 = 4.0);
                     tol = 1e-5, verbose = false)
sol0 = run_to_steady(xuan_shen_config(; Nz = 96, La_t = 0.3, closure = :my25, E6 = 0.0);
                     tol = 1e-5, verbose = false)

# Channel is nondimensional u★=H=1 → νt already = KM/(u★ H)
g = sol4.config.grid
zH = g.zc ./ g.H
KM4 = copy(sol4.state.νt_c)
KM0 = copy(sol0.state.νt_c)
peak4 = maximum(KM4)
peak0 = maximum(KM0)
ratio = peak4 / max(peak0, eps())
σ4 = -zH[argmax(KM4)]
σ0 = -zH[argmax(KM0)]

@printf("  max KM/(u★H)  E6=4 / E6=0 = %.4f / %.4f   ratio=%.2f\n", peak4, peak0, ratio)
@printf("  peak depth σ  E6=4 / E6=0 = %.3f / %.3f\n", σ4, σ0)
ok_E6 = ratio > 1.5
@printf("  [%-4s] E6=4 raises KM by >50%%\n", ok_E6 ? "PASS" : "FAIL")

z_ref, (KM_no, KM_e1, KM_e4) = load_km_csv(joinpath(root, "data", "kc04_fig1_KM.csv"); cols = 2:4)
z_les, (KM_sh, KM_lt) = load_km_csv(joinpath(root, "data", "mcwilliams1997_fig3b_KM.csv"); cols = 2:3)

@printf("\n  Digitized KC04 Fig.1 peaks  E6=4 / noLC = %.3f / %.3f  ratio=%.2f\n",
        maximum(KM_e4), maximum(KM_no), maximum(KM_e4) / maximum(KM_no))
@printf("  Digitized McW LES peaks     La0.3 / shear = %.3f / %.3f\n",
        maximum(KM_lt), maximum(KM_sh))

# ---------------------------------------------------------------------------
# B) McWilliams 短跑：仅作反例说明（不参与硬指标、不写入叠画 MY25 列）
# ---------------------------------------------------------------------------
println()
println("="^72)
println("B) McWilliams+Coriolis short run  [NOT an E6 gate — expect weak Ps]")
println("="^72)

zi = 48.0
cfg_m4 = mcwilliams1997_config(; Nz = 48, H = 90.0, u★ = 0.0061, La_t = 0.3,
                               wavelength = 60.0, f = 1e-4, closure = :my25, E6 = 4.0)
cfg_m0 = mcwilliams1997_config(; Nz = 48, H = 90.0, u★ = 0.0061, La_t = 0.3,
                               wavelength = 60.0, f = 1e-4, closure = :my25, E6 = 0.0)
solm4 = run_to_steady(cfg_m4; tol = 5e-4, max_steps = 8_000, check_every = 200,
                      cfl = 0.15, verbose = false)
solm0 = run_to_steady(cfg_m0; tol = 5e-4, max_steps = 8_000, check_every = 200,
                      cfl = 0.15, verbose = false)
scale = cfg_m4.forcing.u★ * zi
mask(z) = (z .>= -1.0) .& (z .<= 0.0)
zzi4 = solm4.config.grid.zc ./ zi
zzi0 = solm0.config.grid.zc ./ zi
KMm4 = solm4.state.νt_c ./ scale
KMm0 = solm0.state.νt_c ./ scale
m4 = mask(zzi4)
m0 = mask(zzi0)
peakm4 = maximum(KMm4[m4])
peakm0 = maximum(KMm0[m0])
ratio_m = peakm4 / max(peakm0, eps())
@printf("  McW short-run KM/(u★zi) E6=4/0 = %.4f / %.4f  ratio=%.2f\n", peakm4, peakm0, ratio_m)
@printf("  [%-4s] (informational) ratio≈1 means Ps weak — do NOT plot these as E6 proof\n",
        ratio_m < 1.2 ? "WARN" : "NOTE")

# ---------------------------------------------------------------------------
# 写出 CSV：叠画 MY25 = 通道（不是 McW）
# ---------------------------------------------------------------------------
out_csv = joinpath(outdir, "my25_vs_kc04_fig1.csv")
open(out_csv, "w") do io
    println(io, "z_over_zi,KM_KC04_noLC,KM_KC04_E6_1,KM_KC04_E6_4,KM_McWLES_shear,KM_McWLES_La0.3,KM_MY25_channel_E6_0,KM_MY25_channel_E6_4")
    println(io, "# MY25 columns = channel KM/(u★H) on z/H. KC04/McW = KM/(u★ zi). Compare E6 trend, not absolute identity.")
    println(io, "# Do NOT use old column names KM_MY25_E6_* from previous script revisions (those were weak McW runs).")
    for zq in range(0.0, -1.0; length = 81)
        @printf(io, "%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n",
                zq,
                interp1(z_ref, KM_no, zq),
                interp1(z_ref, KM_e1, zq),
                interp1(z_ref, KM_e4, zq),
                interp1(z_les, KM_sh, zq),
                interp1(z_les, KM_lt, zq),
                interp1(zH, KM0, zq),
                interp1(zH, KM4, zq))
    end
end
@printf("\nWrote %s\n", out_csv)

out2 = joinpath(outdir, "my25_channel_E6_compare.csv")
open(out2, "w") do io
    println(io, "z_over_H,nu_t_E6_0,nu_t_E6_4")
    for i in eachindex(g.zc)
        @printf(io, "%.6e,%.6e,%.6e\n", zH[i], KM0[i], KM4[i])
    end
end
@printf("Wrote %s\n", out2)

summary_path = joinpath(outdir, "E6_gate_summary.txt")
open(summary_path, "w") do io
    println(io, "Ocean1DRANS MY25/KC04 E6 gate")
    @printf(io, "channel_peak_E6_4=%.6f\n", peak4)
    @printf(io, "channel_peak_E6_0=%.6f\n", peak0)
    @printf(io, "channel_ratio=%.4f\n", ratio)
    @printf(io, "channel_gate=%s\n", ok_E6 ? "PASS" : "FAIL")
    @printf(io, "mcw_short_ratio=%.4f  (informational; often ~1)\n", ratio_m)
    println(io, "plot_primary=output/my25_channel_E6_compare.csv  columns nu_t_E6_4 vs nu_t_E6_0")
end
@printf("Wrote %s\n", summary_path)

write_profiles_csv(joinpath(outdir, "xuan_shen_La0.3_my25_E6_4.csv"), sol4)
write_profiles_csv(joinpath(outdir, "xuan_shen_La0.3_my25_E6_0.csv"), sol0)

println()
println("="^72)
println("How to plot (E6 proof)")
println("="^72)
println("""
  ★ 主图：output/my25_channel_E6_compare.csv
      x = nu_t_E6_0（细）与 nu_t_E6_4（粗），y = z_over_H
      本跑峰值比 ratio ≈ $(round(ratio; digits = 2))  （须 > 1.5）

  叠画论文：output/my25_vs_kc04_fig1.csv
      KM_KC04_E6_4 / KM_KC04_noLC — 论文 Fig.1
      KM_MY25_channel_E6_*       — 本模型（通道）
      若你仍在画旧列名 KM_MY25_E6_*，那是旧 McW 弱 Ps 结果，会看起来 E6≈无效。
""")

@printf("SUMMARY: %s\n", ok_E6 ? "PASS — channel E6=4 clearly above E6=0" :
                                "FAIL — channel E6 did not enhance KM")
exit(ok_E6 ? 0 : 1)
