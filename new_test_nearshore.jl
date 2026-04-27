# Auto-generated style script for nearshore Langmuir resolvent analysis
# Bottom boundary: solid wall / no-slip
# Surface boundary: rigid-lid stress-free by default

# --- Cell 1 ---
using LinearAlgebra
using Printf
using Statistics
using CairoMakie

# --- Cell 2 ---
const N_demo = 256
const H_demo = 1.0
const USTAR_DEMO = 1.0
const La_t_demo = 0.2
const k0H_demo = 3.5
const NU_T_MEAN_DEMO = 0.07
const NU_WATER_DIM = 1.0e-6

const BOTTOM_BC_DEMO = :no_slip
const SURFACE_BC_DEMO = :stress_free

const U_SURFACE_DEMO = 0.0
const CURRENT_DELTA_OVER_H_DEMO = 0.08
const STOKES_BOTTOM_DAMPING_POWER_DEMO = 1.0

const RECT_COLLOC_STAMP = "2026-04-27-nearshore-solid-wall"

# --- Cell 3: rectangular collocation definitions ---
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

    # Integral weights for Chebyshev-Gauss points mapped from [-1, 1] to [-H, 0].
    w_xi = [(pi / N) * sin((j - 0.5) * pi / N) for j in 1:N]
    w_y_int = (H / 2) .* w_xi

    return RectGrid(N, Nv, Nw, Float64(H), xi_v, xi_w, xi_int, y_v, y_w, y_int,
        Dv, D2v, Dw, D2w, Pv, Pw, I_vw, I_wv, I_Nv_N, I_Nw_N, w_y_int)
end

# --- Cell 4: resolvent operator definitions ---
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

    F12 = (-im * kz) .* (diagc(dUsv) * g.I_vw)
    F21 = (-im * kz) .* (diagc(dUw) * g.I_wv)
    return (; Δv, L_OS, L_Sq, F12, F21)
end

function vbc(coeff_v::AbstractVector, Nv::Int, Nw::Int)
    return vcat(complex.(coeff_v), zeros(ComplexF64, Nw))
end

function wbc(coeff_w::AbstractVector, Nv::Int, Nw::Int)
    return vcat(zeros(ComplexF64, Nv), complex.(coeff_w))
end

function impose_boundary_rows!(M::AbstractMatrix, g::RectGrid, bottom_bc::Symbol, surface_bc::Symbol)
    Nv, Nw, N = g.Nv, g.Nw, g.N
    Dv = complex.(g.Dv)
    D2v = complex.(g.D2v)
    Dw = complex.(g.Dw)

    e_v_bottom = zeros(ComplexF64, Nv); e_v_bottom[1] = 1
    e_v_surface = zeros(ComplexF64, Nv); e_v_surface[end] = 1
    e_w_bottom = zeros(ComplexF64, Nw); e_w_bottom[1] = 1
    e_w_surface = zeros(ComplexF64, Nw); e_w_surface[end] = 1

    # v boundary rows. Bottom y=-H; surface y=0.
    M[N+1, :] .= vbc(e_v_bottom, Nv, Nw)
    if bottom_bc == :no_slip
        M[N+2, :] .= vbc(Dv[1, :], Nv, Nw)
    elseif bottom_bc == :stress_free
        M[N+2, :] .= vbc(D2v[1, :], Nv, Nw)
    else
        error("unknown bottom_bc = $bottom_bc")
    end

    if surface_bc == :no_slip
        M[N+3, :] .= vbc(Dv[end, :], Nv, Nw)
    elseif surface_bc == :stress_free
        M[N+3, :] .= vbc(D2v[end, :], Nv, Nw)
    else
        error("unknown surface_bc = $surface_bc")
    end
    M[N+4, :] .= vbc(e_v_surface, Nv, Nw)

    # omega_y boundary rows.
    if bottom_bc == :no_slip
        M[Nv+N+1, :] .= wbc(e_w_bottom, Nv, Nw)
    elseif bottom_bc == :stress_free
        M[Nv+N+1, :] .= wbc(Dw[1, :], Nv, Nw)
    end

    if surface_bc == :no_slip
        M[Nv+Nw, :] .= wbc(e_w_surface, Nv, Nw)
    elseif surface_bc == :stress_free
        M[Nv+Nw, :] .= wbc(Dw[end, :], Nv, Nw)
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

    impose_boundary_rows!(M, g, bottom_bc, surface_bc)

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

# --- Cell 5: utilities ---
energy_norm_sq(f::AbstractVector, w_y::AbstractVector) = sum(w_y .* abs2.(f))
extract_v_from_uhat(u1::AbstractVector, N::Int) = abs.(u1[N + 1:2N])

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

# --- Cell 6: grid ---
workdir = @__DIR__
g_rect = build_rect_grid(N_demo, H_demo)
@printf("[%s] sum(w_y)=%.12f (should be close to H=%.6f)\n",
    RECT_COLLOC_STAMP, sum(g_rect.w_y_int), H_demo)

# --- Cell 7: nearshore background profiles ---
y_over_H_v = g_rect.y_v ./ H_demo
y_over_H_w = g_rect.y_w ./ H_demo

s_v = @. clamp((g_rect.y_v + H_demo) / H_demo, 0.0, 1.0)
s_w = @. clamp((g_rect.y_w + H_demo) / H_demo, 0.0, 1.0)

# Eulerian mean current: U(-H)=0, U(0)=U_SURFACE_DEMO.
delta_c = max(CURRENT_DELTA_OVER_H_DEMO, 1e-8)
Uv = @. U_SURFACE_DEMO * (1.0 - exp(-s_v / delta_c)) / (1.0 - exp(-1.0 / delta_c))
Uw = @. U_SURFACE_DEMO * (1.0 - exp(-s_w / delta_c)) / (1.0 - exp(-1.0 / delta_c))

# Finite-depth Stokes drift, optionally damped to zero at the solid bottom.
Us_surface = USTAR_DEMO / max(La_t_demo^2, 1e-30)
kh = max(k0H_demo, 1e-8)
Usv = @. Us_surface * cosh(2.0 * kh * s_v) / cosh(2.0 * kh)
Usw = @. Us_surface * cosh(2.0 * kh * s_w) / cosh(2.0 * kh)

p_stokes = max(STOKES_BOTTOM_DAMPING_POWER_DEMO, 0.0)
if p_stokes > 0
    Usv .*= s_v .^ p_stokes
    Usw .*= s_w .^ p_stokes
end

# Eddy viscosity: molecular floor plus parabolic turbulent contribution.
nu_floor = NU_WATER_DIM / max(USTAR_DEMO * H_demo, 1e-30)
nu_mean = max(NU_T_MEAN_DEMO, 0.0)
shape_v = @. max(4.0 * s_v * (1.0 - s_v), 0.0)
shape_w = @. max(4.0 * s_w * (1.0 - s_w), 0.0)
scale_v = nu_mean / max(mean(shape_v), 1e-30)
scale_w = nu_mean / max(mean(shape_w), 1e-30)
nuTv = @. nu_floor + scale_v * shape_v
nuTw = @. nu_floor + scale_w * shape_w

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

ibot_v = argmin(g_rect.y_v)
isurf_v = argmax(g_rect.y_v)

@printf("Profiles: nearshore analytic profiles on Chebyshev grids\n")
@printf("Boundary condition: bottom = %s, surface = %s\n", string(BOTTOM_BC_DEMO), string(SURFACE_BC_DEMO))
@printf("Eulerian U: U(bottom)=%.4e, U(surface)=%.4e\n", Uv[ibot_v], Uv[isurf_v])
@printf("Stokes: finite-depth, La_t=%.3f, k0H=%.4f, bottom damping power=%.2f\n",
    La_t_demo, k0H_demo, p_stokes)
@printf("Stokes values: Us(bottom)=%.4e, Us(surface)=%.4e\n", Usv[ibot_v], Usv[isurf_v])
@printf("nu_T: parabolic turbulent part + molecular floor, target mean = %.4e, nu_floor = %.4e\n",
    nu_mean, nu_floor)
@printf("UL values: UL(bottom)=%.4e, UL(surface)=%.4e\n", ULv[ibot_v], ULv[isurf_v])

# --- Cell 8: background profile plot ---
yH_v = g_rect.y_v ./ H_demo
yH_w = g_rect.y_w ./ H_demo

fig_prof = Figure(size = (520, 900), fontsize = 12)

ax1 = Axis(fig_prof[1, 1]; xlabel = "U^L", ylabel = "y/H", title = "Lagrangian mean U^L")
lines!(ax1, ULv, yH_v; label = "U^L (v)", color = :steelblue, linewidth = 2)
lines!(ax1, ULw, yH_w; label = "U^L (w)", color = :coral, linewidth = 2, linestyle = :dash)
axislegend(ax1; position = :rt)

ax2 = Axis(fig_prof[2, 1]; xlabel = "U^s", ylabel = "y/H", title = "Stokes drift U^s")
lines!(ax2, Usv, yH_v; label = "U^s (v)", color = :darkgreen, linewidth = 2)
lines!(ax2, Usw, yH_w; label = "U^s (w)", color = :olive, linewidth = 2, linestyle = :dash)
axislegend(ax2; position = :rt)

ax3 = Axis(fig_prof[3, 1]; xlabel = "nu_T", ylabel = "y/H", title = "Eddy viscosity nu_T")
lines!(ax3, nuTv, yH_v; label = "nu_T (v)", color = :purple, linewidth = 2)
lines!(ax3, nuTw, yH_w; label = "nu_T (w)", color = :magenta, linewidth = 2, linestyle = :dash)
axislegend(ax3; position = :rt)

fig_prof

# --- Cell 9: solve first resolvent mode ---
kx = 0.0
kz = 2 * pi / H_demo
omega = 1e-8

@printf("Nearshore target mode: La_t=%.1f, (kx*H, kz*H, omega) = (%.4f, %.6f, %.4e)\n",
    La_t_demo, kx * H_demo, kz * H_demo, omega)

res = transfer_gain_rect(
    kx, kz, omega, g_rect,
    Uv, Usv, nuTv, dUv, d2Uv, dUsv, dnuTv, d2nuTv,
    Uw, Usw, nuTw, dUw, d2Uw, dUsw, dnuTw, d2nuTw;
    bottom_bc = BOTTOM_BC_DEMO,
    surface_bc = SURFACE_BC_DEMO,
)

@printf("Euclidean SVD: G = sigma_1^2 = %.6e\n", res.G)

wr = weighted_resolvent_modes(res.T, g_rect.w_y_int)
sigma1 = wr.sigma1
@printf("Weighted SVD: sigma_1 = %.6e, G_w = %.6e\n", sigma1, abs2(sigma1))

psi_phys = sigma1 .* wr.psi1
phi_phys = wr.phi1

psi_vis = copy(psi_phys)
phi_vis = copy(phi_phys)
iy_phase, theta_opt = choose_phase_from_v!(psi_vis, phi_vis, g_rect.N)
@printf("Phase rotation (align v at max|v|, iy=%d): theta = %.4f rad\n", iy_phase, theta_opt)

vprof = extract_v_from_uhat(psi_vis, g_rect.N)
@printf("Check sum w_y|v|^2 (v only) = %.6e\n", energy_norm_sq(vprof, g_rect.w_y_int))

# --- Cell 10: build physical fields in z-y plane ---
py = sortperm(g_rect.y_int)
N = g_rect.N
uh_c = psi_vis[1:N]
vh_c = psi_vis[N + 1:2N]
wh_c = psi_vis[2N + 1:3N]
fx_c = phi_vis[1:N]
fy_c = phi_vis[N + 1:2N]
fz_c = phi_vis[2N + 1:3N]
yh_plot = g_rect.y_int[py] ./ H_demo

nz = 720
zH = range(0.0, 4.0; length = nz)
phys_field(hc::AbstractVector) = [real(hc[i] * cis(kz * (zH[j] * H_demo))) for i in 1:N, j in 1:nz]

Umat = phys_field(uh_c)[py, :]
Vmat = phys_field(vh_c)[py, :]
Wmat = phys_field(wh_c)[py, :]
Dx_mat = phys_field(fx_c)[py, :]
Dy_mat = phys_field(fy_c)[py, :]
Dz_mat = phys_field(fz_c)[py, :]

# --- Cell 11: diagnostic energy profiles ---
w_y_plot = g_rect.w_y_int[py]

uh_raw = psi_vis[1:N][py]
vh_raw = psi_vis[N + 1:2N][py]
wh_raw = psi_vis[2N + 1:3N][py]
fx_raw = phi_vis[1:N][py]
fy_raw = phi_vis[N + 1:2N][py]
fz_raw = phi_vis[2N + 1:3N][py]

w2_raw = abs2.(wh_raw)
dz2_raw = abs2.(fz_raw)
dz2_weighted = w_y_plot .* dz2_raw

fig_diag = Figure(size = (980, 420), fontsize = 12)
ax1 = Axis(fig_diag[1, 1]; xlabel = "value", ylabel = "y/H", title = "Response |w|^2")
lines!(ax1, w2_raw, yh_plot; color = :steelblue, label = "|w|^2 raw", linewidth = 2)
lines!(ax1, w_y_plot .* w2_raw, yh_plot; color = :darkorange, label = "w_y |w|^2", linewidth = 2, linestyle = :dash)
axislegend(ax1; position = :rb)

ax2 = Axis(fig_diag[1, 2]; xlabel = "value", ylabel = "y/H", title = "Forcing |d_z|^2")
lines!(ax2, dz2_raw, yh_plot; color = :seagreen, label = "|d_z|^2 raw", linewidth = 2)
lines!(ax2, dz2_weighted, yh_plot; color = :purple, label = "w_y |d_z|^2", linewidth = 2, linestyle = :dash)
axislegend(ax2; position = :rb)

println("Top boundary raw/weighted d_z^2 ratio = ",
    dz2_raw[end] / max(dz2_weighted[end], 1e-30))
println("Bottom boundary raw/weighted d_z^2 ratio = ",
    dz2_raw[1] / max(dz2_weighted[1], 1e-30))

fig_diag

# --- Cell 12: w and d_z two-dimensional fields ---
fig_wdz = Figure(size = (1100, 450), fontsize = 13)

ax_w = Axis(fig_wdz[1, 1]; xlabel = "z/H", ylabel = "y/H", title = "Response: w")
wmax = max(maximum(abs, Wmat), 1e-12)
ct_w = contourf!(ax_w, collect(zH), yh_plot, Matrix(Wmat');
    levels = range(-wmax, wmax; length = 15),
    extendlow = :auto,
    extendhigh = :auto)
Colorbar(fig_wdz[1, 2], ct_w; label = "w")

ax_dz = Axis(fig_wdz[1, 3]; xlabel = "z/H", ylabel = "y/H", title = "Forcing: d_z")
dzmax = max(maximum(abs, Dz_mat), 1e-12)
ct_dz = contourf!(ax_dz, collect(zH), yh_plot, Matrix(Dz_mat');
    levels = range(-dzmax, dzmax; length = 15),
    extendlow = :auto,
    extendhigh = :auto)
Colorbar(fig_wdz[1, 4], ct_dz; label = "d_z")

fig_wdz

# --- Cell 13: full response and forcing mode plot ---
function add_component_panel!(fig, pos, field, ttl)
    vmax = max(maximum(abs, field), 1e-12)
    ax = Axis(fig[pos...]; xlabel = "z/H", ylabel = "y/H", title = ttl)
    ct = contourf!(ax, collect(zH), yh_plot, Matrix(field');
        levels = range(-vmax, vmax; length = 17),
        extendlow = :auto,
        extendhigh = :auto)
    Colorbar(fig[pos[1], pos[2] + 1], ct)
end

fig_mode = Figure(size = (1400, 900), fontsize = 12)
add_component_panel!(fig_mode, (1, 1), Umat, "Response u")
add_component_panel!(fig_mode, (1, 3), Vmat, "Response v")
add_component_panel!(fig_mode, (1, 5), Wmat, "Response w")
add_component_panel!(fig_mode, (2, 1), Dx_mat, "Forcing d_x")
add_component_panel!(fig_mode, (2, 3), Dy_mat, "Forcing d_y")
add_component_panel!(fig_mode, (2, 5), Dz_mat, "Forcing d_z")

fig_prof_mode = Figure(size = (1100, 420), fontsize = 12)
ax_r = Axis(fig_prof_mode[1, 1]; xlabel = "energy", ylabel = "y/H", title = "Response energy profile")
lines!(ax_r, abs2.(uh_raw), yh_plot; color = :steelblue, label = "u^2 raw", linewidth = 1.3, linestyle = :dash)
lines!(ax_r, abs2.(vh_raw), yh_plot; color = :orangered, label = "v^2 raw", linewidth = 1.3, linestyle = :dash)
lines!(ax_r, abs2.(wh_raw), yh_plot; color = :seagreen, label = "w^2 raw", linewidth = 1.3, linestyle = :dash)
lines!(ax_r, w_y_plot .* abs2.(uh_raw), yh_plot; color = :steelblue, label = "w_y u^2", linewidth = 2)
lines!(ax_r, w_y_plot .* abs2.(vh_raw), yh_plot; color = :orangered, label = "w_y v^2", linewidth = 2)
lines!(ax_r, w_y_plot .* abs2.(wh_raw), yh_plot; color = :seagreen, label = "w_y w^2", linewidth = 2)
axislegend(ax_r; position = :rb)

ax_f = Axis(fig_prof_mode[1, 2]; xlabel = "energy", ylabel = "y/H", title = "Forcing energy profile")
lines!(ax_f, abs2.(fx_raw), yh_plot; color = :steelblue, label = "d_x^2 raw", linewidth = 1.3, linestyle = :dash)
lines!(ax_f, abs2.(fy_raw), yh_plot; color = :orangered, label = "d_y^2 raw", linewidth = 1.3, linestyle = :dash)
lines!(ax_f, abs2.(fz_raw), yh_plot; color = :seagreen, label = "d_z^2 raw", linewidth = 1.3, linestyle = :dash)
lines!(ax_f, w_y_plot .* abs2.(fx_raw), yh_plot; color = :steelblue, label = "w_y d_x^2", linewidth = 2)
lines!(ax_f, w_y_plot .* abs2.(fy_raw), yh_plot; color = :orangered, label = "w_y d_y^2", linewidth = 2)
lines!(ax_f, w_y_plot .* abs2.(fz_raw), yh_plot; color = :seagreen, label = "w_y d_z^2", linewidth = 2)
axislegend(ax_f; position = :rb)

println("sigma1 (weighted) = ", sigma1)

fig_mode
fig_prof_mode
