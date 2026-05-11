# Reproduce the workflow of figure 7 in
# "Resolvent model-based analyses of coherent structures in Langmuir turbulence".
#
# Figure 7 plots the leading-mode energy ratio
#
#     chi_1(kx, kz; c) = sigma_1^2 / sum_j sigma_j^2
#
# for selected phase speeds c = U^L(y_c).  The paper uses depths
# y_c/H = -0.05, -0.12, -0.4, and shows cases La_t = 0.2 and 0.3.
#
# This script computes the model contours.  The grey LES energy-spectrum contours in
# the paper require LES spectral data and are therefore left as an optional overlay.

using LinearAlgebra
using Printf
using Statistics
using CairoMakie
using Serialization
using DelimitedFiles

include(joinpath(@__DIR__, "nearshore_langmuir_resolvent.jl"))

const FIG7_LA_VALUES = (0.2, 0.3)
const FIG7_DEPTHS_OVER_H = (-0.05, -0.12, -0.4)
const FIG7_LAMBDA_X_RANGE = (0.1, 40.0)
const FIG7_LAMBDA_Z_RANGE = (0.05, 40.0)

logrange(a::Real, b::Real, n::Int) = collect(10 .^ range(log10(a), log10(b); length = n))

function weighted_singular_values(
    T::AbstractMatrix,
    w_y::AbstractVector,
)
    wi = 1.0 ./ sqrt.(w_y)
    ws = sqrt.(w_y)
    Winv = Diagonal(vcat(wi, wi, wi))
    Wsqrt = Diagonal(vcat(ws, ws, ws))
    return svdvals(Wsqrt * T * Winv)
end

function leading_energy_ratio_from_transfer(T::AbstractMatrix, w_y::AbstractVector)
    s = weighted_singular_values(T, w_y)
    s2 = abs2.(s)
    return s2[1] / max(sum(s2), 1e-300)
end

function leading_energy_ratio_rect(
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
    return leading_energy_ratio_from_transfer(T, g.w_y_int)
end

function interpolate_profile_at_y(y_grid::AbstractVector, f_grid::AbstractVector, yq::Real)
    p = sortperm(y_grid)
    ys = y_grid[p]
    fs = f_grid[p]
    y = clamp(Float64(yq), ys[1], ys[end])
    j = searchsortedlast(ys, y)
    j = clamp(j, 1, length(ys) - 1)
    t = (y - ys[j]) / (ys[j + 1] - ys[j] + 1e-300)
    return (1 - t) * fs[j] + t * fs[j + 1]
end

function lagrangian_phase_speeds_at_depths(g::RectGrid, profiles, depths_over_H)
    Uv, Usv = profiles[1], profiles[2]
    ULv = Uv .+ Usv
    return [interpolate_profile_at_y(g.y_v, ULv, d * g.H) for d in depths_over_H]
end

function make_fig7_profiles(g::RectGrid, La_t::Real)
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

function write_ratio_csv(path::AbstractString, lambda_x, lambda_z, depths_over_H, ratio)
    rows = Matrix{Float64}(undef, length(lambda_x) * length(lambda_z) * length(depths_over_H), 4)
    r = 1
    for id in eachindex(depths_over_H), ix in eachindex(lambda_x), iz in eachindex(lambda_z)
        rows[r, 1] = depths_over_H[id]
        rows[r, 2] = lambda_x[ix]
        rows[r, 3] = lambda_z[iz]
        rows[r, 4] = ratio[id, ix, iz]
        r += 1
    end
    writedlm(path, rows, ',')
end

function run_fig7_scan(;
    N::Int = 64,
    n_lambda_x::Int = 24,
    n_lambda_z::Int = 28,
    La_values = FIG7_LA_VALUES,
    depths_over_H = FIG7_DEPTHS_OVER_H,
    H::Real = 1.0,
    bottom_bc::Symbol = :stress_free,
    surface_bc::Symbol = :stress_free,
    output_prefix::AbstractString = "fig7_leading_ratio",
    full_paper_grid::Bool = false,
    save_every_lambda_x::Bool = true,
)
    if full_paper_grid
        # The paper does not explicitly restate the Fig.7 sampling count in the caption.
        # Using the same wavelength grid as Fig.3 is a conservative reproduction choice.
        n_lambda_x = 96
        n_lambda_z = 112
    end

    lambda_x = logrange(FIG7_LAMBDA_X_RANGE..., n_lambda_x)
    lambda_z = logrange(FIG7_LAMBDA_Z_RANGE..., n_lambda_z)
    kx_vals = @. 2pi / (lambda_x * H)
    kz_vals = @. 2pi / (lambda_z * H)

    g = build_rect_grid(N, H)
    results = Dict{Float64, Any}()

    @printf("Fig.7 scan setup: N=%d, n_lambda_x=%d, n_lambda_z=%d, n_depth=%d\n",
        N, n_lambda_x, n_lambda_z, length(depths_over_H))
    @printf("Total transfer evaluations per La_t = %d\n",
        n_lambda_x * n_lambda_z * length(depths_over_H))
    @printf("Boundary conditions: bottom=%s, surface=%s\n", string(bottom_bc), string(surface_bc))

    for La_t in La_values
        profiles = make_fig7_profiles(g, La_t)
        c_values = lagrangian_phase_speeds_at_depths(g, profiles, depths_over_H)
        ratio = fill(NaN, length(depths_over_H), n_lambda_x, n_lambda_z)

        @printf("\nLa_t = %.3f\n", La_t)
        for (id, yc) in enumerate(depths_over_H)
            @printf("  depth y/H = %.3f, c = U^L(y) = %.6e\n", yc, c_values[id])
        end

        for ix in eachindex(lambda_x)
            kx = kx_vals[ix]
            for iz in eachindex(lambda_z)
                kz = kz_vals[iz]
                for id in eachindex(depths_over_H)
                    omega = c_values[id] * kx
                    ratio[id, ix, iz] = leading_energy_ratio_rect(
                        kx, kz, omega, g, profiles;
                        bottom_bc = bottom_bc,
                        surface_bc = surface_bc,
                    )
                end
            end

            @printf("  La_t=%.3f progress: ix=%d/%d, lambda_x/H=%.4g\n",
                La_t, ix, length(lambda_x), lambda_x[ix])

            if save_every_lambda_x
                cache_path = "$(output_prefix)_La$(replace(string(La_t), "." => "p"))_cache.jls"
                serialize(cache_path, (; lambda_x, lambda_z, depths_over_H, c_values, ratio, La_t, N))
            end
        end

        csv_path = "$(output_prefix)_La$(replace(string(La_t), "." => "p")).csv"
        write_ratio_csv(csv_path, lambda_x, lambda_z, depths_over_H, ratio)
        @printf("Saved CSV: %s\n", csv_path)

        results[Float64(La_t)] = (; lambda_x, lambda_z, depths_over_H, c_values, ratio)
    end

    return (; g, results, N, n_lambda_x, n_lambda_z, depths_over_H, bottom_bc, surface_bc)
end

function plot_fig7_scan(
    scan;
    output_png::AbstractString = "fig7_leading_mode_ratio.png",
    low_rank_contour::Real = 0.65,
)
    La_values = sort(collect(keys(scan.results)))
    nrow = length(scan.depths_over_H)
    ncol = length(La_values)
    fig = Figure(size = (440 * ncol, 360 * nrow), fontsize = 14)
    levels = range(0.2, 0.95; length = 16)

    for (icol, La_t) in enumerate(La_values)
        r = scan.results[La_t]
        for id in eachindex(scan.depths_over_H)
            ax = Axis(
                fig[id, icol];
                xscale = log10,
                yscale = log10,
                xlabel = id == nrow ? "lambda_x/H" : "",
                ylabel = icol == 1 ? "lambda_z/H" : "",
                title = "La_t=$(La_t), c=U^L(y/H=$(round(scan.depths_over_H[id], digits=2)))",
            )

            cf = contourf!(
                ax,
                r.lambda_x,
                r.lambda_z,
                r.ratio[id, :, :];
                levels = levels,
                colormap = :viridis,
                extendlow = :auto,
                extendhigh = :auto,
            )
            contour!(
                ax,
                r.lambda_x,
                r.lambda_z,
                r.ratio[id, :, :];
                levels = [Float64(low_rank_contour)],
                color = :black,
                linewidth = 2,
            )
            if id == 1
                Colorbar(fig[id, ncol + 1], cf; label = "sigma_1^2 / sum sigma_j^2")
            end
        end
    end

    save(output_png, fig)
    @printf("Saved figure: %s\n", output_png)
    return fig
end

function estimate_fig7_work(;
    N::Int = 256,
    n_lambda_x::Int = 96,
    n_lambda_z::Int = 112,
    n_depth::Int = 3,
    n_La::Int = 2,
)
    n_eval = n_La * n_depth * n_lambda_x * n_lambda_z
    n_state = 2N + 6
    n_inout = 3N
    bytes_T = 16 * n_inout * n_inout
    return (;
        n_eval,
        state_matrix_size = (n_state, n_state),
        transfer_matrix_size = (n_inout, n_inout),
        transfer_matrix_memory_MiB = bytes_T / 2.0^20,
        comment = "Each evaluation performs one dense complex linear solve with $(n_inout) RHS plus a full SVD of a $(n_inout)x$(n_inout) weighted transfer matrix.",
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    # Small default run.  For a paper-scale wavelength grid, call:
    #
    #   scan = run_fig7_scan(N=256, full_paper_grid=true)
    #
    scan = run_fig7_scan()
    plot_fig7_scan(scan)
    println("Paper-sized work estimate: ", estimate_fig7_work())
end
