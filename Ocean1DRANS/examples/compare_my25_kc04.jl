#!/usr/bin/env julia
#=
对照 Kantha & Clayson (2004) Fig.1 / McWilliams et al. (1997) Fig.3b
验证 MY2.5 + KC04 Langmuir 项是否正确实现。

参考数据（手工数字化）：
  data/kc04_fig1_KM.csv
  data/mcwilliams1997_fig3b_KM.csv

运行：
  julia --project=. examples/compare_my25_kc04.jl
=#

using Ocean1DRANS
using Printf

root = joinpath(@__DIR__, "..")
outdir = joinpath(root, "output")
mkpath(outdir)

function load_km_csv(path; zcol=1, cols=2:4)
    zs = Float64[]; cols_data = [Float64[] for _ in cols]
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
        x1, x2 = x[i], x[i+1]
        if (x1 >= xq >= x2) || (x1 <= xq <= x2)
            t = (xq - x1) / (x2 - x1)
            return f[i] + t * (f[i+1] - f[i])
        end
    end
    return f[end]
end

# ---------------------------------------------------------------------------
# 1) 通道算例：应力平衡，干净检验 E6 敏感性（实现正确性的硬指标）
# ---------------------------------------------------------------------------
println("="^72)
println("A) Xuan–Shen channel: MY25 E6=4 vs E6=0 (implementation check)")
println("="^72)

sol4 = run_to_steady(xuan_shen_config(; Nz = 96, La_t = 0.3, closure = :my25, E6 = 4.0);
                     tol = 1e-5, verbose = false)
sol0 = run_to_steady(xuan_shen_config(; Nz = 96, La_t = 0.3, closure = :my25, E6 = 0.0);
                     tol = 1e-5, verbose = false)
peak4 = maximum(sol4.state.νt_c)
peak0 = maximum(sol0.state.νt_c)
ratio = peak4 / max(peak0, eps())
σ4 = -sol4.config.grid.zc[argmax(sol4.state.νt_c)] / sol4.config.grid.H
σ0 = -sol0.config.grid.zc[argmax(sol0.state.νt_c)] / sol0.config.grid.H
@printf("  max νt E6=4/0 = %.4f / %.4f  ratio=%.2f\n", peak4, peak0, ratio)
@printf("  peak depth σ  E6=4/0 = %.3f / %.3f\n", σ4, σ0)
ok_E6 = ratio > 1.5
@printf("  [%-4s] E6=4 raises KM by >50%% (KC04 Fig.1 qualitative)\n", ok_E6 ? "PASS" : "FAIL")

# ---------------------------------------------------------------------------
# 2) McWilliams 1997 设定：与数字化 KC04 Fig.1 / LES Fig.3b 叠画
# ---------------------------------------------------------------------------
println()
println("="^72)
println("B) McWilliams 1997 forcing: MY25 vs digitized KC04 Fig.1")
println("="^72)

# Use mixed-layer scale zi ≈ 48 m (typical McWilliams/KC04 snapshot scale)
zi = 48.0
cfg4 = mcwilliams1997_config(; Nz = 72, H = 90.0, u★ = 0.0061, La_t = 0.3,
                             wavelength = 60.0, f = 1e-4, closure = :my25, E6 = 4.0)
cfg0 = mcwilliams1997_config(; Nz = 72, H = 90.0, u★ = 0.0061, La_t = 0.3,
                             wavelength = 60.0, f = 1e-4, closure = :my25, E6 = 0.0)

# Time-march toward quasi-steady Ekman+Stokes state
solm4 = run_to_steady(cfg4; tol = 5e-4, max_steps = 20_000, check_every = 200,
                      cfl = 0.15, verbose = false)
solm0 = run_to_steady(cfg0; tol = 5e-4, max_steps = 20_000, check_every = 200,
                      cfl = 0.15, verbose = false)

u★ = cfg4.forcing.u★
scale = u★ * zi

function nondim_profile(sol, zi, scale)
    g = sol.config.grid
    zzi = g.zc ./ zi
    KM = sol.state.νt_c ./ scale
    k = sol.state.k ./ u★^2
    return zzi, KM, k
end

zzi4, KM4, k4 = nondim_profile(solm4, zi, scale)
zzi0, KM0, k0 = nondim_profile(solm0, zi, scale)

# Mask to mixed layer z/zi ∈ [-1, 0]
mask4 = (zzi4 .>= -1.0) .& (zzi4 .<= 0.0)
mask0 = (zzi0 .>= -1.0) .& (zzi0 .<= 0.0)
peakm4 = maximum(KM4[mask4])
peakm0 = maximum(KM0[mask0])
σm4 = -zzi4[mask4][argmax(KM4[mask4])]
σm0 = -zzi0[mask0][argmax(KM0[mask0])]
ratio_m = peakm4 / max(peakm0, eps())
@printf("  MY25 McW  max KM/(u*zi) E6=4/0 = %.4f / %.4f  ratio=%.2f\n", peakm4, peakm0, ratio_m)
@printf("  peak σ(= -z/zi) E6=4/0 = %.3f / %.3f   (zi=%.0fm assumed)\n", σm4, σm0, zi)
@printf("  converged E6=4/0: %s / %s  (iters %d / %d)\n",
        solm4.converged, solm0.converged, solm4.iterations, solm0.iterations)

z_ref, (KM_no, KM_e1, KM_e4) = load_km_csv(joinpath(root, "data", "kc04_fig1_KM.csv"); cols = 2:4)
z_les, (KM_sh, KM_lt) = load_km_csv(joinpath(root, "data", "mcwilliams1997_fig3b_KM.csv"); cols = 2:3)

ref_peak4 = maximum(KM_e4)
ref_peak0 = maximum(KM_no)
@printf("  KC04 Fig.1 digitized peaks E6=4 / noLC = %.3f / %.3f  ratio=%.2f\n",
        ref_peak4, ref_peak0, ref_peak4 / ref_peak0)
@printf("  McW LES Fig.3b Langmuir / shear peaks = %.3f / %.3f\n",
        maximum(KM_lt), maximum(KM_sh))

# Order-of-magnitude checks against KC04 Fig.1
ok_mag = 0.08 < peakm4 < 0.45   # E6=4 should be O(0.2), not O(0.03)
ok_enh = peakm4 > peakm0 * 1.2
status_mag = ok_mag ? "PASS" : "WARN"
status_enh = ok_enh ? "PASS" : "WARN"
@printf("  [%-4s] MY25 E6=4 KM/(u*zi) in LES/KC04 order (0.08–0.45)\n", status_mag)
@printf("  [%-4s] MY25 E6=4 KM > E6=0 on McWilliams forcing\n", status_enh)
if !(ok_mag && ok_enh)
    println("  note: McWilliams run may be under-spun (Coriolis + no stratification);")
    println("        treat panel B as overlay reference. Gate on panel A for implementation.")
end

# ---------------------------------------------------------------------------
# 写出叠画用 CSV
# ---------------------------------------------------------------------------
out_csv = joinpath(outdir, "my25_vs_kc04_fig1.csv")
open(out_csv, "w") do io
    println(io, "z_over_zi,KM_KC04_noLC,KM_KC04_E6_1,KM_KC04_E6_4,KM_McWLES_shear,KM_McWLES_La0.3,KM_MY25_E6_0,KM_MY25_E6_4")
    for zq in range(0.0, -1.0; length = 81)
        @printf(io, "%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n",
                zq,
                interp1(z_ref, KM_no, zq),
                interp1(z_ref, KM_e1, zq),
                interp1(z_ref, KM_e4, zq),
                interp1(z_les, KM_sh, zq),
                interp1(z_les, KM_lt, zq),
                interp1(zzi0, KM0, zq),
                interp1(zzi4, KM4, zq))
    end
end
@printf("\nWrote %s\n", out_csv)

# Channel overlay for local E6 test
out2 = joinpath(outdir, "my25_channel_E6_compare.csv")
open(out2, "w") do io
    println(io, "z_over_H,nu_t_E6_0,nu_t_E6_4")
    g = sol4.config.grid
    for i in eachindex(g.zc)
        @printf(io, "%.6e,%.6e,%.6e\n",
                g.zc[i] / g.H, sol0.state.νt_c[i], sol4.state.νt_c[i])
    end
end
@printf("Wrote %s\n", out2)

# Also dump raw McWilliams profiles
write_profiles_csv(joinpath(outdir, "mcwilliams1997_my25_E6_4.csv"), solm4)
write_profiles_csv(joinpath(outdir, "mcwilliams1997_my25_E6_0.csv"), solm0)

println()
println("="^72)
println("How to plot / interpret")
println("="^72)
println("""
  1. 实现是否正确：看 A) — E6=4 必须明显抬高 KM（本仓库要求 ratio>1.5）。
  2. 与论文金标准：画 output/my25_vs_kc04_fig1.csv
       - KM_KC04_E6_4 / KM_KC04_noLC ：Kantha & Clayson (2004) Fig.1
       - KM_McWLES_La0.3 / KM_McWLES_shear ：McWilliams et al. (1997) Fig.3b
       - KM_MY25_E6_4 / KM_MY25_E6_0 ：本模型
  3. 注意：本 McWilliams 配置尚未含浮力分层，zi 取 48 m 作无量纲；
     深部截断不会与 Fig.1 完全一致，应重点比峰值量级与 E6 增强趋势。
  4. 不要用 Xuan–Shen Fig.2 判断 MY25 shape。
""")

passed = ok_E6
@printf("SUMMARY: %s (channel E6 gate); McWilliams overlay written for plotting\n",
        passed ? "MY25/KC04 implementation check PASSED" : "MY25/KC04 implementation check FAILED")
exit(passed ? 0 : 1)
