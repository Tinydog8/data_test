#=
将 Ocean1DRANS 稳态廓线导出为 resolvent 分析可用的格式。

输出列：y/H, U/u★, Us/u★, nu_t/(u★ H)
（与 langmuir_resolvent_cpu.jl 中 LES 数字化剖面约定一致：y∈[0,-H]）
=#

using Ocean1DRANS
using Printf

function export_resolvent_profiles(path::AbstractString, sol::SteadySolution)
    cfg = sol.config
    g = cfg.grid
    u★ = cfg.forcing.u★
    H = g.H
    open(path, "w") do io
        println(io, "y_over_H,U_over_ustar,Us_over_ustar,nu_t_over_ustar_H")
        for i in eachindex(g.zc)
            yH = g.zc[i] / H
            @printf(io, "%.8e,%.8e,%.8e,%.8e\n",
                    yH,
                    sol.state.U[i] / u★,
                    cfg.stokes.us_c[i] / u★,
                    sol.state.νt_c[i] / (u★ * H))
        end
    end
    return path
end

outdir = joinpath(@__DIR__, "..", "output")
mkpath(outdir)

for (La_t, clos, tag) in (
        (0.2, KLStokesClosure(; E6 = 4.0, channel = true), "klstokes"),
        (0.3, KLStokesClosure(; E6 = 4.0, channel = true), "klstokes"),
        (0.2, KPPLTClosure(; Cw = 0.55, use_langmuir = true), "kpplt"),
        (0.3, KPPLTClosure(; Cw = 0.55, use_langmuir = true), "kpplt"),
    )
    cfg = xuan_shen_config(; Nz = 128, La_t = La_t, Reτ = 1000, k0H = 3.5, closure = clos)
    sol = run_to_steady(cfg; tol = 1e-7, verbose = false)
    path = joinpath(outdir, @sprintf("resolvent_base_La%.1f_%s.csv", La_t, tag))
    export_resolvent_profiles(path, sol)
    @printf("La_t=%.1f %-8s  max(νt/(u★H))=%.4f  -> %s\n",
            La_t, tag, maximum(sol.state.νt_c) / (cfg.forcing.u★ * cfg.grid.H), path)
end
