using LinearAlgebra
using Printf
using Statistics

const NEARSHORE_RECT_COLLOC_STAMP = "2026-04-27-nearshore-solid-wall"

function trefethen_cheb_dmat(N::Int)
    N < 2 && error("N >= 2 is required")
    x = [cos(pi * Float64(j) / (N - 1)) for j in 0:N-1]
    c = [i == 1 || i == N ? 2.0 : 1.0 for i in 1:N] .* [(-1)^(i - 1) for i in 1:N]
    D = zeros(Float64, N, N)
    for i in 1:N, j in 1:N
        i == j && continue
        D[i, j] = (c[i] / c[j]) / (x[i] - x[j])
    end
    for i in 1:N
        D[i, i] = -sum(D[i, k] for k in 1:N if k != i)
    end
    return x, D, D * D
end

cheb1_points(N::Int) = [-cos(pi * (j + 0.5) / N) for j in 0:N-1]

function barycentric_weights(n::Int, source_kind::Symbol)
    if source_kind == :lobatto
        w = Float64[((j == 1 || j == n) ? 0.5 : 1.0) * (-1)^(j - 1) for j in 1:n]
    elseif source_kind == :cheb1
        w = Float64[(-1)^(j - 1) * sin((j - 0.5) * pi / n) for j in 1:n]
    else
        error("unknown source_kind = $source_kind")
    end
    w ./= maximum(abs.(w))
    return w
end

function barycentric_interp_matrix(xs::Vector{Float64}, zt::Vector{Float64}; source_kind::Symbol = :lobatto)
    n, m = length(xs), length(zt)
    w = barycentric_weights(n, source_kind)
    P = zeros(ComplexF64, m, n)
    for i in 1:m
        zi = zt[i]
        hit = false
        for j in 1:n
            if abs(zi - xs[j]) < 1e-14
                P[i, :] .= 0
                P[i, j] = 1.0
                hit = true
                break
            end
        end
        hit && continue

        denom = sum(w[k] / (zi - xs[k]) for k in 1:n)
        for j in 1:n
            P[i, j] = (w[j] / (zi - xs[j])) / denom
        end
    end
    return P
end

diagc(v::AbstractVector) = Matrix(Diagonal(complex.(collect(v))))
Delta_hat(D2::AbstractMatrix, k::Real) = complex.(D2) - k^2 * I

mutable struct RectGrid
    N::Int
    Nv::Int
    Nw::Int
    H::Float64
    xi_v::Vector{Float64}
    xi_w::Vector{Float64}
    xi_int::Vector{Float64}
    y_v::Vector{Float64}
    y_w::Vector{Float64}
    y_int::Vector{Float64}
    Dv::Matrix{Float64}
    D2v::Matrix{Float64}
    Dw::Matrix{Float64}
    D2w::Matrix{Float64}
    Pv::Matrix{ComplexF64}
    Pw::Matrix{ComplexF64}
    I_vw::Matrix{ComplexF64}
    I_wv::Matrix{ComplexF64}
    I_Nv_N::Matrix{ComplexF64}
    I_Nw_N::Matrix{ComplexF64}
    w_y_int::Vector{Float64}
end

function build_rect_grid(N_pde::Int, H::Real)
    N_pde >= 8 || error("N_pde >= 8 is required")
    N = N_pde
    Nv, Nw = N + 4, N + 2

    xi_v, Dv, D2v = trefethen_cheb_dmat(Nv)
    xi_w, Dw, D2w = trefethen_cheb_dmat(Nw)
    xi_int = cheb1_points(N)

    sc = -2 / H
    Dv .*= sc
    D2v .*= sc^2
    Dw .*= sc
    D2w .*= sc^2

    y_v = -H .* (xi_v .+ 1) ./ 2
    y_w = -H .* (xi_w .+ 1) ./ 2
    y_int = -H .* (xi_int .+ 1) ./ 2

    Pv = barycentric_interp_matrix(collect(xi_v), collect(xi_int); source_kind = :lobatto)
    Pw = barycentric_interp_matrix(collect(xi_w), collect(xi_int); source_kind = :lobatto)
    I_vw = barycentric_interp_matrix(collect(xi_w), collect(xi_v); source_kind = :lobatto)
    I_wv = barycentric_interp_matrix(collect(xi_v), collect(xi_w); source_kind = :lobatto)
    I_Nv_N = barycentric_interp_matrix(collect(xi_int), collect(xi_v); source_kind = :cheb1)
    I_Nw_N = barycentric_interp_matrix(collect(xi_int), collect(xi_w); source_kind = :cheb1)

    # Chebyshev-Gauss quadrature on x = -cos(theta), mapped from [-1, 1] to [-H, 0].
    w_xi = [(pi / N) * sin((j - 0.5) * pi / N) for j in 1:N]
    w_y_int = (H / 2) .* w_xi

    return RectGrid(N, Nv, Nw, Float64(H), xi_v, xi_w, xi_int, y_v, y_w, y_int,
        Dv, D2v, Dw, D2w, Pv, Pw, I_vw, I_wv, I_Nv_N, I_Nw_N, w_y_int)
end

"""
    make_nearshore_profiles(g; ...)

Construct smooth nearshore background profiles on the v- and omega-grids.

Coordinate convention:
  * y = -H is the bottom wall.
  * y = 0 is the sea surface.

Default physical choices:
  * Eulerian current U(y) satisfies bottom no-slip, U(-H)=0.
  * Stokes drift is a finite-depth analytic profile, optionally damped to zero at the bottom wall.
  * Eddy viscosity has a parabolic turbulent part plus a molecular floor.

The returned tuple is ordered for `transfer_gain_rect`.
"""
function make_nearshore_profiles(
    g::RectGrid;
    ustar::Real = 1.0,
    La_t::Real = 0.2,
    k0H::Real = 3.5,
    U_surface::Real = 1.0,
    U_center::Real = 0.5,
    bottom_roughness_over_H::Real = 1.0e-4,
    surface_roughness_over_H::Real = 1.0e-4,
    log_stitch_center_over_H::Real = 0.5,
    log_blend_half_width_over_H::Real = 0.04,
    nuT_mean::Real = 0.07,
    nu_molecular::Real = 1.0e-6,
    finite_depth_stokes::Bool = true,
    stokes_bottom_damping_power::Real = 1.0,
)
    H = g.H
    k0 = k0H / H

    function wall_coordinate(y)
        return clamp((y + H) / H, 0.0, 1.0)
    end

    function eulerian_current(y)
        s = wall_coordinate(y)
        z0b = max(Float64(bottom_roughness_over_H), 1e-10)
        z0s = max(Float64(surface_roughness_over_H), 1e-10)
        s_match = clamp(Float64(log_stitch_center_over_H), 0.05, 0.95)

        bottom_shape = log((s + z0b) / z0b) / log((s_match + z0b) / z0b)
        surface_shape = log(((1 - s) + z0s) / z0s) / log(((1 - s_match) + z0s) / z0s)

        U_bottom = Float64(U_center) * bottom_shape
        U_surface_branch = Float64(U_surface) - (Float64(U_surface) - Float64(U_center)) * surface_shape

        hw = max(Float64(log_blend_half_width_over_H), 1e-6)
        blend_to_surface = 0.5 * (1 + tanh((s - s_match) / hw))
        return (1 - blend_to_surface) * U_bottom + blend_to_surface * U_surface_branch
    end

    Us0 = Float64(ustar) / max(Float64(La_t)^2, 1e-30)
    function stokes_drift(y)
        s = wall_coordinate(y)
        if finite_depth_stokes
            kh = max(k0H, 1e-8)
            Us = Us0 * cosh(2 * kh * s) / cosh(2 * kh)
        else
            Us = Us0 * exp(2 * k0 * y)
        end
        damping_power = max(Float64(stokes_bottom_damping_power), 0.0)
        damping = damping_power == 0.0 ? 1.0 : s^damping_power
        return Us * damping
    end

    function eddy_viscosity(y)
        s = wall_coordinate(y)
        shape = max(4 * s * (1 - s), 0.0)
        return Float64(nu_molecular) + Float64(nuT_mean) * 1.5 * shape
    end

    Uv = eulerian_current.(g.y_v)
    Uw = eulerian_current.(g.y_w)
    Usv = stokes_drift.(g.y_v)
    Usw = stokes_drift.(g.y_w)
    nuTv = eddy_viscosity.(g.y_v)
    nuTw = eddy_viscosity.(g.y_w)

    dUv = real.(g.Dv * Uv)
    d2Uv = real.(g.D2v * Uv)
    dUsv = real.(g.Dv * Usv)
    dnuTv = real.(g.Dv * nuTv)
    d2nuTv = real.(g.D2v * nuTv)

    dUw = real.(g.Dw * Uw)
    d2Uw = real.(g.D2w * Uw)
    dUsw = real.(g.Dw * Usw)
    dnuTw = real.(g.Dw * nuTw)
    d2nuTw = real.(g.D2w * nuTw)

    return (Uv, Usv, nuTv, dUv, d2Uv, dUsv, dnuTv, d2nuTv,
            Uw, Usw, nuTw, dUw, d2Uw, dUsw, dnuTw, d2nuTw)
end

# Uv,Uw are Eulerian mean currents. Usv,Usw are Stokes drift.
# UL = U + Us is formed internally to avoid double-counting Stokes drift.
function build_L_blocks(kx::Real, kz::Real, g::RectGrid,
    Uv, Usv, nuTv, dUv, d2Uv, dUsv, dnuTv, d2nuTv,
    Uw, Usw, nuTw, dUw, d2Uw, dUsw, dnuTw, d2nuTw)

    Nv, Nw = g.Nv, g.Nw
    k = sqrt(kx^2 + kz^2)
    Δv = Delta_hat(g.D2v, k)
    Δ2v = Δv * Δv
    Δw = Delta_hat(g.D2w, k)

    I_Nv = Matrix{ComplexF64}(I, Nv, Nv)
    I_Nw = Matrix{ComplexF64}(I, Nw, Nw)

    ULv = Uv .+ Usv
    ULw = Uw .+ Usw

    D_ULv = diagc(ULv)
    D_Uppv = diagc(d2Uv)
    D_nuTv = diagc(nuTv)
    D_dnupv = diagc(dnuTv)
    D_d2nupv = diagc(d2nuTv)

    D_ULw = diagc(ULw)
    D_nuTw = diagc(nuTw)
    D_dnupw = diagc(dnuTw)

    Dv_c = complex.(g.Dv)
    D2v_c = complex.(g.D2v)
    Dw_c = complex.(g.Dw)

    L_OS = (-im * kx) .* (D_ULv * Δv) .+ (im * kx) .* D_Uppv .+ (D_nuTv * Δ2v) .+
           (2 .* D_dnupv * Dv_c * Δv) .+ (D_d2nupv * (D2v_c .+ k^2 .* I_Nv))
    L_Sq = (-im * kx) .* D_ULw .+ (D_nuTw * Δw) .+ (D_dnupw * Dw_c)

    # Langmuir vortex-force coupling: Stokes shear enters OS, Eulerian shear enters Squire.
    F12 = (-im * kz) .* (diagc(dUsv) * g.I_vw)
    F21 = (-im * kz) .* (diagc(dUw) * g.I_wv)

    return (; Δv, L_OS, L_Sq, F12, F21)
end

function v_bc_row(coeff_v::AbstractVector, Nv::Int, Nw::Int)
    return vcat(complex.(coeff_v), zeros(ComplexF64, Nw))
end

function omega_bc_row(coeff_w::AbstractVector, Nv::Int, Nw::Int)
    return vcat(zeros(ComplexF64, Nv), complex.(coeff_w))
end

function set_v_boundary_rows!(M::AbstractMatrix, g::RectGrid, bottom_bc::Symbol, surface_bc::Symbol)
    Nv, Nw, N = g.Nv, g.Nw, g.N
    Dv = complex.(g.Dv)
    D2v = complex.(g.D2v)
    e_bottom = zeros(ComplexF64, Nv); e_bottom[1] = 1
    e_surface = zeros(ComplexF64, Nv); e_surface[end] = 1

    # Bottom wall at y = -H.
    M[N + 1, :] .= v_bc_row(e_bottom, Nv, Nw)
    if bottom_bc == :no_slip
        M[N + 2, :] .= v_bc_row(Dv[1, :], Nv, Nw)       # Dv = 0 gives tangential no-slip.
    elseif bottom_bc == :stress_free
        M[N + 2, :] .= v_bc_row(D2v[1, :], Nv, Nw)
    else
        error("unsupported bottom_bc = $bottom_bc")
    end

    # Sea surface at y = 0.
    if surface_bc == :stress_free
        M[N + 3, :] .= v_bc_row(D2v[end, :], Nv, Nw)
        M[N + 4, :] .= v_bc_row(e_surface, Nv, Nw)
    elseif surface_bc == :no_slip
        M[N + 3, :] .= v_bc_row(Dv[end, :], Nv, Nw)
        M[N + 4, :] .= v_bc_row(e_surface, Nv, Nw)
    else
        error("unsupported surface_bc = $surface_bc")
    end
    return M
end

function set_omega_boundary_rows!(M::AbstractMatrix, g::RectGrid, bottom_bc::Symbol, surface_bc::Symbol)
    Nv, Nw, N = g.Nv, g.Nw, g.N
    Dw = complex.(g.Dw)
    e_bottom = zeros(ComplexF64, Nw); e_bottom[1] = 1
    e_surface = zeros(ComplexF64, Nw); e_surface[end] = 1

    if bottom_bc == :no_slip
        M[Nv + N + 1, :] .= omega_bc_row(e_bottom, Nv, Nw) # omega_y = 0 at a no-slip wall.
    elseif bottom_bc == :stress_free
        M[Nv + N + 1, :] .= omega_bc_row(Dw[1, :], Nv, Nw)
    else
        error("unsupported bottom_bc = $bottom_bc")
    end

    if surface_bc == :stress_free
        M[Nv + Nw, :] .= omega_bc_row(Dw[end, :], Nv, Nw)
    elseif surface_bc == :no_slip
        M[Nv + Nw, :] .= omega_bc_row(e_surface, Nv, Nw)
    else
        error("unsupported surface_bc = $surface_bc")
    end
    return M
end

function build_M_B_C_rect(kx::Real, kz::Real, omega::Real, g::RectGrid,
    Uv, Usv, nuTv, dUv, d2Uv, dUsv, dnuTv, d2nuTv,
    Uw, Usw, nuTw, dUw, d2Uw, dUsw, dnuTw, d2nuTw;
    bottom_bc::Symbol = :no_slip,
    surface_bc::Symbol = :stress_free)

    N, Nv, Nw, Pv, Pw = g.N, g.Nv, g.Nw, g.Pv, g.Pw
    k = sqrt(kx^2 + kz^2)
    k == 0 && error("kx and kz cannot both be zero")

    bl = build_L_blocks(kx, kz, g, Uv, Usv, nuTv, dUv, d2Uv, dUsv, dnuTv, d2nuTv,
        Uw, Usw, nuTw, dUw, d2Uw, dUsw, dnuTw, d2nuTw)
    Δv, L_OS, L_Sq, F12, F21 = bl.Δv, bl.L_OS, bl.L_Sq, bl.F12, bl.F21

    I_Nw = Matrix{ComplexF64}(I, Nw, Nw)
    ntot = Nv + Nw
    M = zeros(ComplexF64, ntot, ntot)

    M[1:N, 1:Nv] = Pv * (im * omega .* Δv .- L_OS)
    M[1:N, Nv+1:end] = -Pv * F12
    M[Nv+1:Nv+N, 1:Nv] = -Pw * F21
    M[Nv+1:Nv+N, Nv+1:end] = Pw * (im * omega .* I_Nw .- L_Sq)

    set_v_boundary_rows!(M, g, bottom_bc, surface_bc)
    set_omega_boundary_rows!(M, g, bottom_bc, surface_bc)

    I_Nv = Matrix{ComplexF64}(I, Nv, Nv)
    Znw = zeros(ComplexF64, Nw, Nw)
    Dv = complex.(g.Dv)

    Binner = hcat((-im * kx) .* Dv, (-k^2) .* I_Nv, (-im * kz) .* Dv)
    Bv = Pv * Binner * kron(Matrix{ComplexF64}(I, 3, 3), g.I_Nv_N)
    Bw = Pw * hcat(im * kz .* I_Nw, Znw, (-im * kx) .* I_Nw) *
         kron(Matrix{ComplexF64}(I, 3, 3), g.I_Nw_N)

    Bmat = zeros(ComplexF64, ntot, 3N)
    Bmat[1:N, :] .= Bv
    Bmat[Nv+1:Nv+N, :] .= Bw

    Ivw = g.I_vw
    row_u_v = (im * kx / k^2) .* (Pv * Dv)
    row_u_w = (-im * kz / k^2) .* (Pv * Ivw)
    row_v_v = Pv
    Znvw = zeros(ComplexF64, N, Nw)
    row_w_v = (im * kz / k^2) .* (Pv * Dv)
    row_w_w = (im * kx / k^2) .* (Pv * Ivw)
    Cmat = vcat(hcat(row_u_v, row_u_w),
                hcat(row_v_v, Znvw),
                hcat(row_w_v, row_w_w))

    return M, Bmat, Cmat
end

function transfer_gain_rect(kx::Real, kz::Real, omega::Real, g::RectGrid,
    Uv, Usv, nuTv, dUv, d2Uv, dUsv, dnuTv, d2nuTv,
    Uw, Usw, nuTw, dUw, d2Uw, dUsw, dnuTw, d2nuTw;
    bottom_bc::Symbol = :no_slip,
    surface_bc::Symbol = :stress_free)

    M, Bmat, Cmat = build_M_B_C_rect(kx, kz, omega, g,
        Uv, Usv, nuTv, dUv, d2Uv, dUsv, dnuTv, d2nuTv,
        Uw, Usw, nuTw, dUw, d2Uw, dUsw, dnuTw, d2nuTw;
        bottom_bc = bottom_bc,
        surface_bc = surface_bc)
    X = M \ Bmat
    T = Cmat * X
    s1 = svdvals(T)[1]
    return (; G = abs2(s1), s1, T, M, Cmat, g, bottom_bc, surface_bc)
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

function choose_phase_from_v!(psi::AbstractVector, phi::AbstractVector, N::Int)
    vh = @view psi[N+1:2N]
    iy = argmax(abs.(vh))
    theta = -angle(vh[iy])
    rot = cis(theta)
    psi .*= rot
    phi .*= rot
    return iy, theta
end

function modal_energy_profiles(psi::AbstractVector, phi::AbstractVector, g::RectGrid)
    N = g.N
    py = sortperm(g.y_int)
    yh = g.y_int[py] ./ g.H
    wy = g.w_y_int[py]

    u = psi[1:N][py]
    v = psi[N+1:2N][py]
    w = psi[2N+1:3N][py]
    dx = phi[1:N][py]
    dy = phi[N+1:2N][py]
    dz = phi[2N+1:3N][py]

    return (;
        y_over_H = yh,
        weights = wy,
        response_raw = (; u = abs2.(u), v = abs2.(v), w = abs2.(w)),
        forcing_raw = (; dx = abs2.(dx), dy = abs2.(dy), dz = abs2.(dz)),
        response_weighted = (; u = wy .* abs2.(u), v = wy .* abs2.(v), w = wy .* abs2.(w)),
        forcing_weighted = (; dx = wy .* abs2.(dx), dy = wy .* abs2.(dy), dz = wy .* abs2.(dz)),
    )
end

function run_nearshore_demo(;
    N::Int = 256,
    H::Real = 1.0,
    kx::Real = 0.0,
    kz::Real = 2pi,
    omega::Real = 1e-8,
    ustar::Real = 1.0,
    La_t::Real = 0.2,
    k0H::Real = 3.5,
    U_surface::Real = 1.0,
    U_center::Real = 0.5,
    bottom_roughness_over_H::Real = 1.0e-4,
    surface_roughness_over_H::Real = 1.0e-4,
    log_stitch_center_over_H::Real = 0.5,
    log_blend_half_width_over_H::Real = 0.04,
    nuT_mean::Real = 0.07,
    nu_molecular::Real = 1e-6,
    stokes_bottom_damping_power::Real = 1.0,
)
    g = build_rect_grid(N, H)
    profiles = make_nearshore_profiles(g;
        ustar = ustar,
        La_t = La_t,
        k0H = k0H,
        U_surface = U_surface,
        U_center = U_center,
        bottom_roughness_over_H = bottom_roughness_over_H,
        surface_roughness_over_H = surface_roughness_over_H,
        log_stitch_center_over_H = log_stitch_center_over_H,
        log_blend_half_width_over_H = log_blend_half_width_over_H,
        nuT_mean = nuT_mean,
        nu_molecular = nu_molecular,
        stokes_bottom_damping_power = stokes_bottom_damping_power)

    res = transfer_gain_rect(kx, kz, omega, g, profiles...;
        bottom_bc = :no_slip,
        surface_bc = :stress_free)
    wr = weighted_resolvent_modes(res.T, g.w_y_int)

    psi = wr.sigma1 .* wr.psi1
    phi = copy(wr.phi1)
    iy, theta = choose_phase_from_v!(psi, phi, g.N)
    energy = modal_energy_profiles(psi, phi, g)

    @printf("[%s]\n", NEARSHORE_RECT_COLLOC_STAMP)
    @printf("Boundary conditions: bottom = no_slip, surface = stress_free\n")
    @printf("sum(w_y) = %.12e, H = %.12e\n", sum(g.w_y_int), H)
    @printf("kxH = %.6e, kzH = %.6e, omega = %.6e\n", kx * H, kz * H, omega)
    @printf("Euclidean sigma1 = %.6e, G = %.6e\n", res.s1, res.G)
    @printf("Weighted sigma1 = %.6e, G_w = %.6e\n", wr.sigma1, abs2(wr.sigma1))
    @printf("Phase aligned at v index %d, theta = %.6e rad\n", iy, theta)
    @printf("max raw |dz|^2 = %.6e\n", maximum(energy.forcing_raw.dz))
    @printf("max weighted w_y|dz|^2 = %.6e\n", maximum(energy.forcing_weighted.dz))

    return (; g, profiles, res, wr, psi, phi, energy)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_nearshore_demo()
end
