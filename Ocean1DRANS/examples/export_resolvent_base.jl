#!/usr/bin/env julia
#=
导出 resolvent 可用的基流（与 Xuan & Shen 2025 一致）：
  UL ≈ Us（解析指数）,  νt 来自 LES 数字化 Fig.2b
=#

using Ocean1DRANS
using Printf

function export_resolvent_profiles(path::AbstractString, sol::SteadySolution)
    cfg = sol.config
    g = cfg.grid
    u★ = cfg.forcing.u★
    H = g.H
    open(path, "w") do io
        println(io, "y_over_H,U_over_ustar,UL_over_ustar,Us_over_ustar,nu_t_over_ustar_H")
        for i in eachindex(g.zc)
            yH = g.zc[i] / H
            UL = sol.state.U[i] + cfg.stokes.us_c[i]
            @printf(io, "%.8e,%.8e,%.8e,%.8e,%.8e\n",
                    yH,
                    sol.state.U[i] / u★,
                    UL / u★,
                    cfg.stokes.us_c[i] / u★,
                    sol.state.νt_c[i] / (u★ * H))
        end
    end
    return path
end

outdir = joinpath(@__DIR__, "..", "output")
mkpath(outdir)

for La_t in (0.2, 0.3)
    cfg = xuan_shen_config(; Nz = 128, La_t = La_t, Reτ = 1000, k0H = 3.5, closure = :les)
    sol = run_to_steady(cfg; verbose = false)
    path = joinpath(outdir, @sprintf("resolvent_base_La%.1f_les.csv", La_t))
    export_resolvent_profiles(path, sol)
    @printf("La_t=%.1f  max|U|/max|Us|=%.3f  max(νt/(u★H))=%.4f  -> %s\n",
            La_t,
            maximum(abs, sol.state.U) / maximum(abs, cfg.stokes.us_c),
            maximum(sol.state.νt_c) / (cfg.forcing.u★ * cfg.grid.H),
            path)
end
