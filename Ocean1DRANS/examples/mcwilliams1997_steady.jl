#!/usr/bin/env julia
#=
McWilliams et al. (1997) 型开洋 Langmuir 混合层稳态算例（含 Coriolis）。
=#

using Ocean1DRANS
using Printf

outdir = joinpath(@__DIR__, "..", "output")
mkpath(outdir)

cfg = mcwilliams1997_config(;
    Nz = 64,
    H = 90.0,
    u★ = 0.0061,
    La_t = 0.3,
    wavelength = 60.0,
    f = 1e-4,
    E6 = 4.0,
)

sol = run_to_steady(cfg; tol = 5e-5, max_steps = 400_000, check_every = 1000, verbose = true)
csv = joinpath(outdir, "mcwilliams1997_La0.3.csv")
write_profiles_csv(csv, sol)
@printf("Wrote %s\n", csv)
@printf("surface U,V = (%.4f, %.4f) m/s, max νt = %.4e m²/s\n",
        sol.state.U[end], sol.state.V[end], maximum(sol.state.νt_c))
