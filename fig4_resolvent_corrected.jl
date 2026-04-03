using LinearAlgebra
using Printf
using DelimitedFiles
using CairoMakie

include(joinpath(@__DIR__, "rectangular_collocation.jl"))

# Paper-consistent notes:
# 1) The resolvent operator is built around the turbulent mean state from LES.
#    For Fig. 4, do not replace the mean flow with U = 0 unless you are making
#    an explicit idealised test unrelated to the paper.
# 2) The paper uses the LES-derived eddy viscosity profile directly; there is no
#    near-surface manual taper such as forcing nu_T -> 0 above a chosen z-cut.
# 3) The energy-optimal modes must be extracted from the weighted SVD associated
#    with the discrete version of the inner product in (2.21), not from an
#    ordinary Euclidean SVD of T.

const N_demo = 256
const H_demo = 1.0
const USTAR_DEMO = 1.0
const La_t_demo = 0.2
const k0H_demo = 3.5
const VELOCITY_CSV_IS_LAGRANGIAN = true

energy_norm_sq(f::AbstractVector, w_y::AbstractVector) = sum(w_y .* abs2.(f))
extract_v_from_uhat(u1::AbstractVector, N::Int) = abs.(u1[N + 1:2N])

function read_csv_pair_resolvent(path::AbstractString)
    m = readdlm(path, ',', Float64)
    size(m, 2) >= 2 || error("需要至少两列: $path")
    return vec(m[:, 1]), vec(m[:, 2])
end

function sort_merge_mean_resolvent(z::Vector{Float64}, v::Vector{Float64})
    length(z) == length(v) || error("z 与 v 长度不一致")
    p = sortperm(z)
    zs = z[p]
    vs = v[p]
    outz = Float64[]
    outv = Float64[]
    i = 1
    n = length(zs)
    eps_z = 1e-12 * max(1.0, maximum(zs) - minimum(zs))
    while i <= n
        z0 = zs[i]
        s = vs[i]
        c = 1
        j = i + 1
        while j <= n && abs(zs[j] - z0) <= eps_z
            s += vs[j]
            c += 1
            j += 1
        end
        push!(outz, z0)
        push!(outv, s / c)
        i = j
    end
    return outz, outv
end

function interp_clamp_resolvent(xsrc::Vector{Float64}, vsrc::Vector{Float64}, xdst::AbstractVector)
    length(xsrc) >= 2 || error("插值点过少")
    p = sortperm(xsrc)
    xs = xsrc[p]
    vs = vsrc[p]
    out = similar(xdst, Float64)
    for k in eachindex(xdst)
        x = xdst[k]
        if x <= xs[1]
            out[k] = vs[1]
        elseif x >= xs[end]
            out[k] = vs[end]
        else
            j = searchsortedlast(xs, x)
            x1, x2 = xs[j], xs[j + 1]
            a = (x - x1) / (x2 - x1)
            out[k] = (1 - a) * vs[j] + a * vs[j + 1]
        end
    end
    return out
end

function resample_uniform_y(y_sorted::AbstractVector, f_sorted::AbstractVector, n::Int = 480)
    yq = range(y_sorted[1], y_sorted[end]; length = n)
    out = zeros(Float64, n)
    for k in eachindex(yq)
        y = yq[k]
        j = searchsortedlast(y_sorted, y)
        j = clamp(j, 1, length(y_sorted) - 1)
        y1, y2 = y_sorted[j], y_sorted[j + 1]
        t = (y - y1) / (y2 - y1 + 1e-30)
        out[k] = (1 - t) * f_sorted[j] + t * f_sorted[j + 1]
    end
    return collect(yq), out
end

function load_profile_on_grid(csv_path::AbstractString, ydst::AbstractVector)
    raw_a, raw_b = read_csv_pair_resolvent(csv_path)
    z_src, val_src = sort_merge_mean_resolvent(raw_a, raw_b)
    return interp_clamp_resolvent(z_src, val_src, ydst), z_src, val_src
end

function chebyshev_polynomial_matrix(y::AbstractVector, degree::Int, H::Real)
    x = @. 2.0 * y / H + 1.0
    M = zeros(Float64, length(y), degree + 1)
    M[:, 1] .= 1.0
    if degree >= 1
        M[:, 2] .= x
    end
    for n in 2:degree
        @views M[:, n + 1] .= 2.0 .* x .* M[:, n] .- M[:, n - 1]
    end
    return M
end

function smooth_profile_chebfit(
    z_src::AbstractVector,
    v_src::AbstractVector,
    ydst::AbstractVector,
    H::Real;
    degree::Int = 24,
)
    deg = min(degree, max(length(z_src) - 1, 0))
    deg >= 1 || error("用于拟合的源点太少")
    A = chebyshev_polynomial_matrix(z_src, deg, H)
    c = A \ Float64.(v_src)
    B = chebyshev_polynomial_matrix(ydst, deg, H)
    return B * c
end

function weighted_resolvent_modes(T::AbstractMatrix, w_y::AbstractVector)
    wi = 1.0 ./ sqrt.(w_y)
    ws = sqrt.(w_y)
    Winv = Diagonal(vcat(wi, wi, wi))
    Wsqrt = Diagonal(vcat(ws, ws, ws))
    Ttw = Wsqrt * T * Winv
    Sw = svd(Ttw)
    sigma1 = Sw.S[1]
    phi1 = Winv * Sw.V[:, 1]
    psi1 = Winv * Sw.U[:, 1]
    return (; sigma1, phi1, psi1, Sw)
end

function choose_phase_from_v!(psi_vis::AbstractVector, phi_vis::AbstractVector, N::Int)
    vh = @view psi_vis[N + 1:2N]
    iy_phase = argmax(abs.(vh))
    theta_opt = -angle(vh[iy_phase])
    rot = cis(theta_opt)
    psi_vis .*= rot
    phi_vis .*= rot
    return iy_phase, theta_opt
end

function main()
    workdir = @__DIR__
    vel_csv = joinpath(workdir, "vel_resolvent.csv")
    nu_csv = joinpath(workdir, "viscosity_resolvent.csv")
    isfile(vel_csv) || error("缺少 vel_resolvent.csv。要复现论文 Fig.4，必须使用 LES 提取的欧拉平均流 U(y)，不能把 U 设为 0。")
    isfile(nu_csv) || error("缺少 viscosity_resolvent.csv")

    g_rect = build_rect_grid(N_demo, H_demo)
    @printf("[%s] sum(w_y)=%.12f (应≈H=%.6f)\n", RECT_COLLOC_STAMP, sum(g_rect.w_y_int), H_demo)

    _, z_u, UL_src = load_profile_on_grid(vel_csv, g_rect.y_v)
    _, z_nu, nu_src = load_profile_on_grid(nu_csv, g_rect.y_v)

    # Paper §2.2: prescribed deep-water monochromatic Stokes drift.
    y_over_H_v = g_rect.y_v ./ H_demo
    y_over_H_w = g_rect.y_w ./ H_demo
    Usv = USTAR_DEMO .* ((1 / La_t_demo^2) .* exp.(2 * k0H_demo .* y_over_H_v))
    Usw = USTAR_DEMO .* ((1 / La_t_demo^2) .* exp.(2 * k0H_demo .* y_over_H_w))

    # Paper §2.3 reports and plots the Lagrangian mean velocity U^L.
    # The linear operator still needs both U^L for advection and the Eulerian
    # mean U (through U' and U''). Reconstruct U from U^L - U^s by default.
    if VELOCITY_CSV_IS_LAGRANGIAN
        Uv_from_csv = smooth_profile_chebfit(z_u, UL_src, g_rect.y_v, H_demo) .- Usv
        Uw_from_csv = smooth_profile_chebfit(z_u, UL_src, g_rect.y_w, H_demo) .- Usw
    else
        Uv_from_csv = smooth_profile_chebfit(z_u, UL_src, g_rect.y_v, H_demo)
        Uw_from_csv = smooth_profile_chebfit(z_u, UL_src, g_rect.y_w, H_demo)
    end

    nuTv = smooth_profile_chebfit(z_nu, nu_src, g_rect.y_v, H_demo)
    nuTw = smooth_profile_chebfit(z_nu, nu_src, g_rect.y_w, H_demo)

    Uv = copy(Uv_from_csv)
    Uw = copy(Uw_from_csv)
    ULv = Uv .+ Usv
    ULw = Uw .+ Usw

    dUv = g_rect.Dv * Uv
    d2Uv = g_rect.D2v * Uv
    dUsv = g_rect.Dv * Usv
    dnuTv = g_rect.Dv * nuTv
    d2nuTv = g_rect.D2v * nuTv

    dUw = g_rect.Dw * Uw
    d2Uw = g_rect.D2w * Uw
    dUsw = g_rect.Dw * Usw
    dnuTw = g_rect.Dw * nuTw
    d2nuTw = g_rect.D2w * nuTw

    isurf_v = argmax(g_rect.y_v)
    @printf("Loaded %d velocity points, %d nu_T points\n", length(z_u), length(z_nu))
    @printf("Velocity CSV interpreted as %s\n", VELOCITY_CSV_IS_LAGRANGIAN ? "U^L(y)" : "Eulerian U(y)")
    @printf("Profiles are Chebyshev-smoothed before differentiation to suppress boundary noise in U'', nu_T'', and derived modes\n")
    @printf("If the CSVs were digitised from plotted curves rather than exported LES arrays, noticeable discrepancies from the paper's Fig.4 can remain even with a correct discretisation.\n")
    @printf("Surface values: U(y=0)=%.4f, U^s(y=0)=%.4f, U^L(y=0)=%.4f\n",
        Uv[isurf_v], Usv[isurf_v], ULv[isurf_v])
    @printf("Derivative diagnostics: max|U'|=%.4e, max|U''|=%.4e, max|nu_T'|=%.4e, max|nu_T''|=%.4e\n",
        maximum(abs, dUv), maximum(abs, d2Uv), maximum(abs, dnuTv), maximum(abs, d2nuTv))

    kx = 0.0
    kz = 2 * pi / H_demo
    omega = 0.0
    @printf("Fig.4 target mode: La_t=%.1f, (kx*H, kz*H, omega) = (%.4f, %.6f, %.4f)\n",
        La_t_demo, kx * H_demo, kz * H_demo, omega)

    res = transfer_gain_rect(
        kx, kz, omega, g_rect,
        Uv, Usv, nuTv, dUv, d2Uv, dUsv, dnuTv, d2nuTv,
        Uw, Usw, nuTw, dUw, d2Uw, dUsw, dnuTw, d2nuTw,
    )
    @printf("Euclidean SVD: G = sigma_1^2 = %.6e\n", res.G)

    wr = weighted_resolvent_modes(res.T, g_rect.w_y_int)
    sigma1 = wr.sigma1
    @printf("Weighted SVD (paper 2.20-2.21): sigma_1 = %.6e, G_w = %.6e\n", sigma1, abs2(sigma1))

    # Paper (2.22): T * phi1 = sigma1 * psi1.
    psi_phys = sigma1 .* wr.psi1
    phi_phys = wr.phi1

    clim_u = 0.8
    clim_d = 2.4
    psi_vis = copy(psi_phys)
    phi_vis = copy(phi_phys)
    psi_vis .*= clim_u / max(maximum(abs.(real.(psi_vis))), 1e-30)
    phi_vis .*= clim_d / max(maximum(abs.(real.(phi_vis))), 1e-30)

    iy_phase, theta_opt = choose_phase_from_v!(psi_vis, phi_vis, g_rect.N)
    @printf("Phase rotation (align Re(v) at max|v|, iy=%d): theta = %.4f rad\n", iy_phase, theta_opt)

    vprof = extract_v_from_uhat(psi_vis, g_rect.N)
    @printf("Check sum w_y|v|^2 (v only) = %.6e\n", energy_norm_sq(vprof, g_rect.w_y_int))

    py = sortperm(g_rect.y_int)
    N = g_rect.N
    uh_c = psi_vis[1:N]
    vh_c = psi_vis[N + 1:2N]
    wh_c = psi_vis[2N + 1:3N]
    fx_c = phi_vis[1:N]
    fy_c = phi_vis[N + 1:2N]
    fz_c = phi_vis[2N + 1:3N]
    yh_plot = g_rect.y_int[py] ./ H_demo

    nz = 768
    zH = range(0.0, 4.0; length = nz)
    phys_field(hc::AbstractVector) = [real(hc[i] * cis(kz * (zH[j] * H_demo))) for i in 1:N, j in 1:nz]

    Umat = phys_field(uh_c)[py, :]
    Vmat = phys_field(vh_c)[py, :]
    Wmat = phys_field(wh_c)[py, :]
    Dx_mat = phys_field(fx_c)[py, :]
    Dy_mat = phys_field(fy_c)[py, :]
    Dz_mat = phys_field(fz_c)[py, :]

    umx = max(maximum(abs, Umat), 1e-30)
    dmx = max(maximum(abs, Dx_mat), 1e-30)
    lev_u = range(-umx, umx; length = 15)
    lev_d = range(-dmx, dmx; length = 15)

    fig4 = Figure(size = (1100, 900), fontsize = 13)

    ax_a = Axis(fig4[1, 1]; xlabel = "z/H", ylabel = "y/H", title = "(a) Response: u, arrows (w,v)")
    ct_a = contourf!(ax_a, collect(zH), yh_plot, Matrix(Umat'); levels = lev_u, extendlow = :auto, extendhigh = :auto)
    Colorbar(fig4[1, 2], ct_a; label = "u")

    Ns = size(Umat, 1)
    step_i = max(1, div(Ns, 18))
    step_j = max(1, div(nz, 28))
    xs_a = Float32[]
    ys_a = Float32[]
    us_a = Float32[]
    vs_a = Float32[]
    len_arrow = Float32(0.045 * H_demo)
    for ii in 1:step_i:Ns
        for j in 1:step_j:nz
            wx, vy = Wmat[ii, j], Vmat[ii, j]
            nm = Float32(hypot(wx, vy) + 1e-20)
            push!(xs_a, Float32(zH[j]))
            push!(ys_a, Float32(yh_plot[ii]))
            push!(us_a, Float32(wx / nm))
            push!(vs_a, Float32(vy / nm))
        end
    end
    arrows2d!(ax_a, xs_a, ys_a, us_a, vs_a;
        lengthscale = len_arrow, normalize = true, color = (:black, 0.45),
        shaftwidth = 1.5f0, tipwidth = 6f0, tiplength = 5f0)

    ax_b = Axis(fig4[1, 3]; xlabel = "|component|^2", ylabel = "y/H", title = "(b) Response |u_i|^2")
    yq, u2q = resample_uniform_y(yh_plot, abs2.(uh_c[py]))
    _, v2q = resample_uniform_y(yh_plot, abs2.(vh_c[py]))
    _, w2q = resample_uniform_y(yh_plot, abs2.(wh_c[py]))
    lines!(ax_b, u2q, yq; label = "u^2", color = :steelblue, linewidth = 2)
    lines!(ax_b, v2q, yq; label = "v^2", color = :orangered, linewidth = 2)
    lines!(ax_b, w2q, yq; label = "w^2", color = :seagreen, linewidth = 2)
    axislegend(ax_b; position = :rt)

    ax_c = Axis(fig4[2, 1]; xlabel = "z/H", ylabel = "y/H", title = "(c) Forcing: d_x, arrows (d_z,d_y)")
    ct_c = contourf!(ax_c, collect(zH), yh_plot, Matrix(Dx_mat'); levels = lev_d, extendlow = :auto, extendhigh = :auto)
    Colorbar(fig4[2, 2], ct_c; label = "d_x")
    xs_c = Float32[]
    ys_c = Float32[]
    us_c = Float32[]
    vs_c = Float32[]
    for ii in 1:step_i:Ns
        for j in 1:step_j:nz
            dzv, dyv = Dz_mat[ii, j], Dy_mat[ii, j]
            nm = Float32(hypot(dzv, dyv) + 1e-20)
            push!(xs_c, Float32(zH[j]))
            push!(ys_c, Float32(yh_plot[ii]))
            push!(us_c, Float32(dzv / nm))
            push!(vs_c, Float32(dyv / nm))
        end
    end
    arrows2d!(ax_c, xs_c, ys_c, us_c, vs_c;
        lengthscale = len_arrow, normalize = true, color = (:black, 0.45),
        shaftwidth = 1.5f0, tipwidth = 6f0, tiplength = 5f0)

    ax_d = Axis(fig4[2, 3]; xlabel = "|component|^2", ylabel = "y/H", title = "(d) Forcing |d_i|^2")
    yqd, dx2q = resample_uniform_y(yh_plot, abs2.(fx_c[py]), 960)
    _, dy2q = resample_uniform_y(yh_plot, abs2.(fy_c[py]), 960)
    _, dz2q = resample_uniform_y(yh_plot, abs2.(fz_c[py]), 960)
    lines!(ax_d, dx2q, yqd; label = "d_x^2", color = :steelblue, linewidth = 2)
    lines!(ax_d, dy2q, yqd; label = "d_y^2", color = :orangered, linewidth = 2)
    lines!(ax_d, dz2q, yqd; label = "d_z^2", color = :seagreen, linewidth = 2)
    axislegend(ax_d; position = :rt)

    save(joinpath(workdir, "fig4_resolvent_corrected.png"), fig4)

    fig_wdz = Figure(size = (1100, 450), fontsize = 13)
    ax_w = Axis(fig_wdz[1, 1]; xlabel = "z/H", ylabel = "y/H", title = "Response: w")
    wmax = max(maximum(abs, Wmat), 1e-12)
    ct_w = contourf!(ax_w, collect(zH), yh_plot, Matrix(Wmat'); levels = range(-wmax, wmax; length = 15), extendlow = :auto, extendhigh = :auto)
    Colorbar(fig_wdz[1, 2], ct_w; label = "w")

    ax_dz = Axis(fig_wdz[1, 3]; xlabel = "z/H", ylabel = "y/H", title = "Forcing: d_z")
    dzmax = max(maximum(abs, Dz_mat), 1e-12)
    ct_dz = contourf!(ax_dz, collect(zH), yh_plot, Matrix(Dz_mat'); levels = range(-dzmax, dzmax; length = 15), extendlow = :auto, extendhigh = :auto)
    Colorbar(fig_wdz[1, 4], ct_dz; label = "d_z")
    save(joinpath(workdir, "fig4_resolvent_w_dz_corrected.png"), fig_wdz)

    return nothing
end

main()
