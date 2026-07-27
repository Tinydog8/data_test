#!/usr/bin/env julia
#=
对照 Kantha & Clayson (2004) Fig.1 / McWilliams et al. (1997) Fig.3b

说明（为何以前 E6≈无效、且低于 KC04）：
  1. 旧脚本用时间推进 + H=90 m + 底拖曳，惯性振荡未收敛，Ps 很弱。
  2. 应用通道结果顶替叠画是错误的——Fig.1 是开洋混合层，不是通道。
  3. 正确设定：zi=33 m、稳态 Ekman–Stokes、仅表面壁面律、预后 q²/q²ℓ。
  4. Kantha et al. (2010) 指出 KC04 原文 E6=4 为笔误，物理值 ≈ 7.2；
     用 E6=7.2 才能对上 Fig.1 粗红线量级（数字化峰值 ≈ 0.22）。

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

function nondim_KM(sol, zi)
    u★ = sol.config.forcing.u★
    return sol.config.grid.zc ./ zi, sol.state.νt_c ./ (u★ * zi)
end

println("="^72)
println("MY25 vs KC04 Fig.1  (McWilliams zi=33 m, steady Ekman + prognostic q²ℓ)")
println("="^72)

zi = 33.0
# KC04 "without Langmuir terms" = Ps≡0 in both turb equations
sol_no = run_to_steady(mcwilliams1997_config(; Nz = 72, H = zi, closure = :my25,
                                              E6 = 0.0, stokes_production = false);
                       tol = 5e-5, max_steps = 400, underrelax = 0.3, verbose = false)
sol_e4 = run_to_steady(mcwilliams1997_config(; Nz = 72, H = zi, closure = :my25, E6 = 4.0);
                       tol = 5e-5, max_steps = 400, underrelax = 0.3, verbose = false)
sol_e72 = run_to_steady(mcwilliams1997_config(; Nz = 72, H = zi, closure = :my25, E6 = 7.2);
                        tol = 1e-4, max_steps = 500, underrelax = 0.3, verbose = false)

zzi_no, KM_no_m = nondim_KM(sol_no, zi)
zzi_e4, KM_e4_m = nondim_KM(sol_e4, zi)
zzi_e72, KM_e72_m = nondim_KM(sol_e72, zi)

peak_no = maximum(KM_no_m)
peak_e4 = maximum(KM_e4_m)
peak_e72 = maximum(KM_e72_m)

@printf("  MY25 noLC     max KM/(u★zi) = %.4f   (KC04 digitized ≈ 0.102)\n", peak_no)
@printf("  MY25 E6=4     max KM/(u★zi) = %.4f   ratio/noLC = %.2f\n", peak_e4, peak_e4 / peak_no)
@printf("  MY25 E6=7.2   max KM/(u★zi) = %.4f   ratio/noLC = %.2f  (Kantha 2010 correction)\n",
        peak_e72, peak_e72 / peak_no)
@printf("  converged no/E6=4/E6=7.2: %s / %s / %s\n",
        sol_no.converged, sol_e4.converged, sol_e72.converged)

z_ref, (KM_no_ref, KM_e1_ref, KM_e4_ref) =
    load_km_csv(joinpath(root, "data", "kc04_fig1_KM.csv"); cols = 2:4)
z_les, (KM_sh, KM_lt) =
    load_km_csv(joinpath(root, "data", "mcwilliams1997_fig3b_KM.csv"); cols = 2:3)

@printf("\n  Digitized KC04 Fig.1  E6=4 / noLC = %.3f / %.3f  ratio=%.2f\n",
        maximum(KM_e4_ref), maximum(KM_no_ref), maximum(KM_e4_ref) / maximum(KM_no_ref))
@printf("  Digitized McW LES     La0.3 / shear = %.3f / %.3f\n",
        maximum(KM_lt), maximum(KM_sh))

ok_no = 0.07 < peak_no < 0.14
ok_e72 = peak_e72 > peak_no * 1.5 && 0.15 < peak_e72 < 0.35
ok_e4_trend = peak_e4 > peak_no * 1.05
@printf("\n  [%-4s] noLC peak in KC04 order (0.07–0.14)\n", ok_no ? "PASS" : "FAIL")
@printf("  [%-4s] E6=4 raises KM vs noLC\n", ok_e4_trend ? "PASS" : "FAIL")
@printf("  [%-4s] E6=7.2 matches Fig.1 thick-line order (≳1.5× noLC, peak 0.15–0.35)\n",
        ok_e72 ? "PASS" : "FAIL")

# Overlay CSV — MY25 from McWilliams/KC04 setup (not channel)
out_csv = joinpath(outdir, "my25_vs_kc04_fig1.csv")
open(out_csv, "w") do io
    println(io, "z_over_zi,KM_KC04_noLC,KM_KC04_E6_1,KM_KC04_E6_4,KM_McWLES_shear,KM_McWLES_La0.3,KM_MY25_noLC,KM_MY25_E6_4,KM_MY25_E6_7p2")
    println(io, "# All KM/(u★ zi). MY25 from McWilliams zi=33m steady Ekman+prognostic q2l.")
    println(io, "# E6=7.2 is Kantha et al. (2010) correction; aligns with digitized KC04 'E6=4' curve.")
    for zq in range(0.0, -1.0; length = 81)
        @printf(io, "%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n",
                zq,
                interp1(z_ref, KM_no_ref, zq),
                interp1(z_ref, KM_e1_ref, zq),
                interp1(z_ref, KM_e4_ref, zq),
                interp1(z_les, KM_sh, zq),
                interp1(z_les, KM_lt, zq),
                interp1(zzi_no, KM_no_m, zq),
                interp1(zzi_e4, KM_e4_m, zq),
                interp1(zzi_e72, KM_e72_m, zq))
    end
end
@printf("\nWrote %s\n", out_csv)

# Also keep channel E6 gate (algebraic, implementation check independent of McW)
println()
println("="^72)
println("Channel E6 gate (algebraic MY25; separate from KC04 Fig.1)")
println("="^72)
solc4 = run_to_steady(xuan_shen_config(; Nz = 64, La_t = 0.3, closure = :my25, E6 = 4.0);
                      tol = 1e-5, verbose = false)
solc0 = run_to_steady(xuan_shen_config(; Nz = 64, La_t = 0.3, closure = :my25, E6 = 0.0);
                      tol = 1e-5, verbose = false)
rc = maximum(solc4.state.νt_c) / max(maximum(solc0.state.νt_c), eps())
@printf("  channel max νt E6=4/0 = %.4f / %.4f  ratio=%.2f\n",
        maximum(solc4.state.νt_c), maximum(solc0.state.νt_c), rc)
ok_ch = rc > 1.5
@printf("  [%-4s] channel E6=4 ≫ E6=0\n", ok_ch ? "PASS" : "FAIL")

out2 = joinpath(outdir, "my25_channel_E6_compare.csv")
open(out2, "w") do io
    println(io, "z_over_H,nu_t_E6_0,nu_t_E6_4")
    g = solc4.config.grid
    for i in eachindex(g.zc)
        @printf(io, "%.6e,%.6e,%.6e\n", g.zc[i] / g.H, solc0.state.νt_c[i], solc4.state.νt_c[i])
    end
end
@printf("Wrote %s\n", out2)

write_profiles_csv(joinpath(outdir, "mcwilliams1997_my25_noLC.csv"), sol_no)
write_profiles_csv(joinpath(outdir, "mcwilliams1997_my25_E6_4.csv"), sol_e4)
write_profiles_csv(joinpath(outdir, "mcwilliams1997_my25_E6_7p2.csv"), sol_e72)

summary_path = joinpath(outdir, "E6_gate_summary.txt")
open(summary_path, "w") do io
    println(io, "Ocean1DRANS MY25 / KC04 Fig.1")
    @printf(io, "mcw_noLC_peak=%.4f\n", peak_no)
    @printf(io, "mcw_E6_4_peak=%.4f\n", peak_e4)
    @printf(io, "mcw_E6_7p2_peak=%.4f\n", peak_e72)
    @printf(io, "kc04_digitized_E6_4=%.4f\n", maximum(KM_e4_ref))
    @printf(io, "channel_ratio_E6_4_over_0=%.4f\n", rc)
    @printf(io, "gate=%s\n", (ok_no && ok_e72 && ok_ch) ? "PASS" : "FAIL")
end

println()
println("="^72)
println("How to plot")
println("="^72)
println("""
  对 KC04 Fig.1：output/my25_vs_kc04_fig1.csv
    谁和谁比：
      KM_MY25_noLC   ↔ KM_KC04_noLC     （无 Langmuir）
      KM_MY25_E6_7p2 ↔ KM_KC04_E6_4     （有 Langmuir；7.2=原文笔误更正）
      KM_MY25_E6_4   仅作对照（原文印刷值，增强偏弱）
    不要拿通道曲线去比 Fig.1。
""")

ok = ok_no && ok_e72 && ok_ch
@printf("SUMMARY: %s\n", ok ? "PASS" : "FAIL")
exit(ok ? 0 : 1)
