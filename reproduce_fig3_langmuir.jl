# Reproduce the workflow of figure 3 in
# "Resolvent model-based analyses of coherent structures in Langmuir turbulence".
#
# The paper computes
#
#     Gmax(kx, kz) = max_omega G(kx, kz, omega)
#
# over wavelength ranges
#
#     lambda_x/H in [0.1, 40], lambda_z/H in [0.05, 40],
#
# using 96 x 112 wavelength samples and 50 phase-speed samples
#
#     c in [0.01 * Umax^L, Umax^L], omega = c * kx.
#
# This script uses the current collocation/resolvent implementation in
# `nearshore_langmuir_resolvent.jl`.  Without the LES mean profiles used in the
# paper, the default background profile is the analytic nearshore profile from
# that file.  The scan parameters below include a small default smoke-test grid;
# set `full_paper_grid = true` or pass larger values to `run_fig3_scan` for the
# paper-sized sweep.

using LinearAlgebra
using Printf
using Statistics
using CairoMakie
using Serialization
using DelimitedFiles

include(joinpath(@__DIR__, "nearshore_langmuir_resolvent.jl"))

const FIG3_LA_VALUES = (0.2, 0.3)
const FIG3_LAMBDA_X_RANGE = (0.1, 40.0)
const FIG3_LAMBDA_Z_RANGE = (0.05, 40.0)
const FIG3_PHASE_SPEED_FRACTION_RANGE = (0.01, 1.0)

logrange(a::Real, b::Real, n::Int) = collect(10 .^ range(log10(a), log10(b); length = n))

function weighted_gain_from_transfer(T::AbstractMatrix, w_y::AbstractVector)
    wi = 1.0 ./ sqrt.(w_y)
    ws = sqrt.(w_y)
    Winv = Diagonal(vcat(wi, wi, wi))
    Wsqrt = Diagonal(vcat(ws, ws, ws))
    sigma1 = svdvals(Wsqrt * T * Winv)[1]
    return abs2(sigma1)
end

function weighted_gain_rect(
    kx::Real,
    kz::Real,
    omega::Real,
    g::RectGrid,
    profiles;
    bottom_bc::Symbol = :stress_free,
    surface_bc::Symbol = :stress_free,
)
    M, Bmat, Cmat = build_M_B_C_rect(
        kx, kz, omega, g, profiles...;
        bottom_bc = bottom_bc,
        surface_bc = surface_bc,
    )
    T = Cmat * (M \ Bmat)
    return weighted_gain_from_transfer(T, g.w_y_int)
end

function lagrangian_mean_max(profiles)
    Uv, Usv = profiles[1], profiles[2]
    Uw, Usw = profiles[9], profiles[10]
    return maximum(vcat(Uv .+ Usv, Uw .+ Usw))
end

function make_fig3_profiles(g::RectGrid, La_t::Real)
    return make_nearshore_profiles(
        g;
        ustar = 1.0,
        La_t = La_t,
        k0H = 3.5,
        U_surface = 1.0,
        U_center = 0.5,
        bottom_roughness_over_H = 1.0e-4,
        surface_roughness_over_H = 1.0e-4,
        log_stitch_center_over_H = 0.5,
        log_blend_half_width_over_H = 0.04,
        nuT_mean = 0.07,
        nu_molecular = 1.0e-6,
        finite_depth_stokes = true,
        stokes_bottom_damping_power = 1.0,
    )
end

function write_scan_csv(path::AbstractString, lambda_x, lambda_z, Gmax, c_at_max)
    rows = Matrix{Float64}(undef, length(lambda_x) * length(lambda_z), 4)
    r = 1
    for ix in eachindex(lambda_x), iz in eachindex(lambda_z)
        rows[r, 1] = lambda_x[ix]
        rows[r, 2] = lambda_z[iz]
        rows[r, 3] = Gmax[ix, iz]
        rows[r, 4] = c_at_max[ix, iz]
        r += 1
    end
    writedlm(path, rows, ',')
end

function run_fig3_scan(;
    N::Int = 64,
    n_lambda_x::Int = 24,
    n_lambda_z::Int = 28,
    n_phase_speed::Int = 12,
    La_values = FIG3_LA_VALUES,
    H::Real = 1.0,
    bottom_bc::Symbol = :stress_free,
    surface_bc::Symbol = :stress_free,
    output_prefix::AbstractString = "fig3_resolvent",
    full_paper_grid::Bool = false,
    save_every_lambda_x::Bool = true,
)
    if full_paper_grid
        n_lambda_x = 96
        n_lambda_z = 112
        n_phase_speed = 50
    end

    lambda_x = logrange(FIG3_LAMBDA_X_RANGE..., n_lambda_x)
    lambda_z = logrange(FIG3_LAMBDA_Z_RANGE..., n_lambda_z)
    kx_vals = @. 2pi / (lambda_x * H)
    kz_vals = @. 2pi / (lambda_z * H)

    g = build_rect_grid(N, H)
    results = Dict{Float64, Any}()

    @printf("Fig.3 scan setup: N=%d, n_lambda_x=%d, n_lambda_z=%d, n_c=%d\n",
        N, n_lambda_x, n_lambda_z, n_phase_speed)
    @printf("Total transfer evaluations per La_t = %d\n",
        n_lambda_x * n_lambda_z * n_phase_speed)
    @printf("Boundary conditions: bottom=%s, surface=%s\n", string(bottom_bc), string(surface_bc))

    for La_t in La_values
        profiles = make_fig3_profiles(g, La_t)
        UmaxL = lagrangian_mean_max(profiles)
        c_values = collect(range(
            FIG3_PHASE_SPEED_FRACTION_RANGE[1] * UmaxL,
            FIG3_PHASE_SPEED_FRACTION_RANGE[2] * UmaxL;
            length = n_phase_speed,
        ))

        Gmax = fill(NaN, n_lambda_x, n_lambda_z)
        c_at_max = fill(NaN, n_lambda_x, n_lambda_z)
        omega_at_max = fill(NaN, n_lambda_x, n_lambda_z)

        @printf("\nLa_t = %.3f: Umax^L = %.6e, c in [%.6e, %.6e]\n",
            La_t, UmaxL, first(c_values), last(c_values))

        for ix in eachindex(lambda_x)
            kx = kx_vals[ix]
            for iz in eachindex(lambda_z)
                kz = kz_vals[iz]
                best_G = -Inf
                best_c = NaN
                best_omega = NaN

                for c in c_values
                    omega = c * kx
                    G = weighted_gain_rect(
                        kx, kz, omega, g, profiles;
                        bottom_bc = bottom_bc,
                        surface_bc = surface_bc,
                    )
                    if G > best_G
                        best_G = G
                        best_c = c
                        best_omega = omega
                    end
                end

                Gmax[ix, iz] = best_G
                c_at_max[ix, iz] = best_c
                omega_at_max[ix, iz] = best_omega
            end

            @printf("  La_t=%.3f progress: ix=%d/%d, lambda_x/H=%.4g\n",
                La_t, ix, length(lambda_x), lambda_x[ix])

            if save_every_lambda_x
                cache_path = "$(output_prefix)_La$(replace(string(La_t), "." => "p"))_cache.jls"
                serialize(cache_path, (; lambda_x, lambda_z, Gmax, c_at_max, omega_at_max, La_t, N))
            end
        end

        csv_path = "$(output_prefix)_La$(replace(string(La_t), "." => "p")).csv"
        write_scan_csv(csv_path, lambda_x, lambda_z, Gmax, c_at_max)
        results[Float64(La_t)] = (; lambda_x, lambda_z, Gmax, c_at_max, omega_at_max, UmaxL)
        @printf("Saved CSV: %s\n", csv_path)
    end

    return (; g, results, N, n_lambda_x, n_lambda_z, n_phase_speed, bottom_bc, surface_bc)
end

function plot_fig3_scan(scan; output_png::AbstractString = "fig3_resolvent_Gmax.png")
    La_values = sort(collect(keys(scan.results)))
    ncols = length(La_values)
    fig = Figure(size = (560 * ncols, 470), fontsize = 14)

    all_logG = Float64[]
    for La_t in La_values
        append!(all_logG, vec(log10.(max.(scan.results[La_t].Gmax, 1e-300))))
    end
    finite_logG = all_logG[isfinite.(all_logG)]
    lo = floor(minimum(finite_logG))
    hi = ceil(maximum(finite_logG))
    levels = range(lo, hi; length = 17)

    for (icol, La_t) in enumerate(La_values)
        r = scan.results[La_t]
        ax = Axis(
            fig[1, icol];
            xscale = log10,
            yscale = log10,
            xlabel = "lambda_x/H",
            ylabel = "lambda_z/H",
            title = "La_t = $(La_t)",
        )
        cf = contourf!(
            ax,
            r.lambda_x,
            r.lambda_z,
            log10.(max.(r.Gmax, 1e-300));
            levels = levels,
            colormap = :viridis,
            extendlow = :auto,
            extendhigh = :auto,
        )
        line_min = max(minimum(r.lambda_x), minimum(r.lambda_z))
        line_max = min(maximum(r.lambda_x), maximum(r.lambda_z))
        lines!(ax, [line_min, line_max], [line_min, line_max]; color = :white, linestyle = :dash, linewidth = 2)
        Colorbar(fig[1, ncols + icol], cf; label = "log10 Gmax")
    end

    save(output_png, fig)
    @printf("Saved figure: %s\n", output_png)
    return fig
end

function estimate_fig3_work(;
    N::Int = 256,
    n_lambda_x::Int = 96,
    n_lambda_z::Int = 112,
    n_phase_speed::Int = 50,
    n_La::Int = 2,
)
    n_eval = n_La * n_lambda_x * n_lambda_z * n_phase_speed
    n_state = 2N + 6
    n_inout = 3N
    bytes_T = 16 * n_inout * n_inout
    return (;
        n_eval,
        state_matrix_size = (n_state, n_state),
        transfer_matrix_size = (n_inout, n_inout),
        transfer_matrix_memory_MiB = bytes_T / 2.0^20,
        comment = "Each evaluation performs one dense complex linear solve with $(n_inout) RHS plus one SVD of a $(n_inout)x$(n_inout) transfer matrix.",
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    # Small default run.  For the paper-sized grid, call:
    #
    #   scan = run_fig3_scan(N=256, full_paper_grid=true)
    #
    scan = run_fig3_scan()
    plot_fig3_scan(scan)
    println("Paper-sized work estimate: ", estimate_fig3_work())
end
