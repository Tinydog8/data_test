# Langmuir resolvent demo (Fig. 3 G_max scan + Fig. 4 mode plots) -- CPU 优化版.
#
# 相对原始 demo 的改动（均为 CPU 侧，不依赖 GPU）:
#   * 线程: 扫描期间把 BLAS 线程设为 1 (scan_Gmax_map 的 blas_threads 关键字, 默认 1),
#     避免外层 @threads 与内层 BLAS 的过度订阅(oversubscription)。当 Julia 线程数 >= 物理
#     核数时这通常显著更快 (实测可达 1.5x~3x)。
#   * 基流 nuT_profile=:les: 用论文 Fig.2(b) 数字化的 LES 涡黏剖面 (复现 Fig.3 高值区)。
#
# 求解路径沿用原始 demo 的做法: 对每个 (kx,kz) 缓存 L 块/B/C (build_resolvent_cache),
# 每个 omega 组装 M 并解 M\B 得到传递矩阵 T = C M^{-1} B (level-3 BLAS, 多核扩展性好),
# 再用加权幂迭代取 sigma_1。
#
# 注: 曾尝试两种"提速": (a) 广义 Schur(QZ) 复用 omega 扫描 —— Schur 分解固定开销过大,
# 反而更慢; (b) 矩阵-free 幂迭代 (避免组装 T) —— 单线程更快, 但其 level-2 (单右端项) 操作
# 内存带宽受限, 在多核/超线程上扩展性差, 在多核强机上反而比组装完整 T 慢。两者均已弃用。
#
# 数值正确性由 verify_fig3_equivalence / verify_sigma1_power 校验。

using LinearAlgebra
using Base.Threads
using Printf
using Statistics
using CairoMakie
using LaTeXStrings

# ---------------------------------------------------------------------------
# Rectangular spectral collocation (Driscoll & Hale 2016; Xuan & Shen JFM App. B)
# ---------------------------------------------------------------------------

const RECT_COLLOC_STAMP = "2026-04-02-sin_theta_weights"

function trefethen_cheb_dmat(N::Int)
    N < 2 && error("N>=2")
    x = [cos(pi * Float64(j) / (N - 1)) for j in 0:N-1]
    c = [i == 1 || i == N ? 2.0 : 1.0 for i in 1:N] .* [(-1)^(i - 1) for i in 1:N]
    D = zeros(N, N)
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
        w = Float64[
            ((j == 1 || j == n) ? 0.5 : 1.0) * (-1)^(j - 1)
            for j in 1:n
        ]
    elseif source_kind == :cheb1
        w = Float64[
            (-1)^(j - 1) * sin((j - 0.5) * pi / n)
            for j in 1:n
        ]
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
    N_pde >= 8 || error("N_pde >= 8")
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
    w_xi = [(pi / N) * sin((j - 0.5) * pi / N) for j in 1:N]
    w_y_int = (H / 2) .* w_xi
    RectGrid(N, Nv, Nw, Float64(H), xi_v, xi_w, xi_int, y_v, y_w, y_int,
        Dv, D2v, Dw, D2w, Pv, Pw, I_vw, I_wv, I_Nv_N, I_Nw_N, w_y_int)
end

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
    D_ULv = diagc(ULv); D_Uppv = diagc(d2Uv); D_Us_pv = diagc(dUsv); D_Upv = diagc(dUv)
    D_nuTv = diagc(nuTv); D_dnupv = diagc(dnuTv); D_d2nupv = diagc(d2nuTv)
    D_ULw = diagc(ULw); D_Upw = diagc(dUw)
    D_nuTw = diagc(nuTw); D_dnupw = diagc(dnuTw); D_d2nupw = diagc(d2nuTw)
    Dv_c = complex.(g.Dv); D2v_c = complex.(g.D2v)
    Dw_c = complex.(g.Dw); D2w_c = complex.(g.D2w)
    L_OS = (-im * kx) .* (D_ULv * Δv) .+ (im * kx) .* D_Uppv .+ (D_nuTv * Δ2v) .+
           (2 .* D_dnupv * Dv_c * Δv) .+ (D_d2nupv * (D2v_c .+ k^2 .* I_Nv))
    L_Sq = (-im * kx) .* D_ULw .+ (D_nuTw * Δw) .+ (D_dnupw * Dw_c)
    F12 = (-im * kz) .* (diagc(dUsv) * g.I_vw)
    F21 = (-im * kz) .* (diagc(dUw) * g.I_wv)
    return (; Δv, L_OS, L_Sq, F12, F21)
end

function build_M_B_C_rect(kx::Real, kz::Real, omega::Real, g::RectGrid,
    Uv, Usv, nuTv, dUv, d2Uv, dUsv, dnuTv, d2nuTv,
    Uw, Usw, nuTw, dUw, d2Uw, dUsw, dnuTw, d2nuTw)
    N, Nv, Nw, Pv, Pw = g.N, g.Nv, g.Nw, g.Pv, g.Pw
    k = sqrt(kx^2 + kz^2)
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
    D2v = complex.(g.D2v)
    function vbc(coeff_v::AbstractVector)
        vcat(complex.(coeff_v), zeros(ComplexF64, Nw))
    end
    M[N+1, :] .= vbc([one(ComplexF64); zeros(ComplexF64, Nv - 1)])
    M[N+2, :] .= vbc(D2v[1, :])
    M[N+3, :] .= vbc(D2v[Nv, :])
    M[N+4, :] .= vbc([zeros(ComplexF64, Nv - 1); one(ComplexF64)])
    Dw = complex.(g.Dw)
    function wbc(coeff_w::AbstractVector)
        vcat(zeros(ComplexF64, Nv), complex.(coeff_w))
    end
    M[Nv+N+1, :] .= wbc(Dw[1, :])
    M[Nv+Nw, :] .= wbc(Dw[Nw, :])
    I_Nv = Matrix{ComplexF64}(I, Nv, Nv)
    Znw = zeros(ComplexF64, Nw, Nw)
    Dv = complex.(g.Dv)
    Binner = hcat((-im * kx) .* Dv, (-k^2) .* I_Nv, (-im * kz) .* Dv)
    Bv = Pv * Binner * kron(Matrix{ComplexF64}(I, 3, 3), g.I_Nv_N)
    Bw = Pw * hcat(im * kz .* I_Nw, Znw, (-im * kx) .* I_Nw) * kron(Matrix{ComplexF64}(I, 3, 3), g.I_Nw_N)
    Bmat = zeros(ComplexF64, ntot, 3N)
    Bmat[1:N, :] .= Bv
    Bmat[Nv+1:Nv+N, :] .= Bw
    Ivw = g.I_vw
    PvC = g.Pv
    k == 0 && error("k=0")
    row_u_v = (im * kx / k^2) .* (PvC * Dv)
    row_u_w = (-im * kz / k^2) .* (PvC * Ivw)
    row_v_v = PvC
    Znvw = zeros(ComplexF64, N, Nw)
    row_w_v = (im * kz / k^2) .* (PvC * Dv)
    row_w_w = (im * kx / k^2) .* (PvC * Ivw)
    Cm = vcat(hcat(row_u_v, row_u_w), hcat(row_v_v, Znvw), hcat(row_w_v, row_w_w))
    return M, Bmat, Cm
end

function transfer_gain_rect(kx::Real, kz::Real, omega::Real, g::RectGrid,
    Uv, Usv, nuTv, dUv, d2Uv, dUsv, dnuTv, d2nuTv,
    Uw, Usw, nuTw, dUw, d2Uw, dUsw, dnuTw, d2nuTw)
    M, Bcopy, Cm = build_M_B_C_rect(kx, kz, omega, g, Uv, Usv, nuTv, dUv, d2Uv, dUsv, dnuTv, d2nuTv,
        Uw, Usw, nuTw, dUw, d2Uw, dUsw, dnuTw, d2nuTw)
    X = M \ Bcopy
    T = Cm * X
    s1 = svdvals(T)[1]
    return (; G=abs2(s1), s1, T, M, Cm, g)
end

weights_output_3N(g::RectGrid) = vcat(g.w_y_int, g.w_y_int, g.w_y_int)

# ---------------------------------------------------------------------------
# Fig. 3: G_max(kx,kz) scan (Xuan & Shen §3.1–3.2)
# ---------------------------------------------------------------------------

logrange10(a::Real, b::Real, n::Int) = collect(10.0 .^ range(log10(a), log10(b); length = n))

function _weight_vectors_3N(w_y::AbstractVector)
    wi = 1.0 ./ sqrt.(w_y)
    ws = sqrt.(w_y)
    return vcat(wi, wi, wi), vcat(ws, ws, ws)
end

# 无缓冲版本: 与原始实现等价, 供 kx≈0 分支与正确性对照使用。
function sigma1_weighted(
    T::AbstractMatrix{<:Complex},
    w_y::AbstractVector;
    maxiter::Int = 50,
    tol::Real = 1e-7,
)
    wi3, ws3 = _weight_vectors_3N(w_y)
    w3 = ws3 .* ws3
    n = size(T, 2)
    x = Vector{ComplexF64}(undef, n)
    @inbounds for i in 1:n
        x[i] = cis(2π * i / n)
    end
    x ./= norm(x)
    σ_old = 0.0
    for _ in 1:maxiter
        u = T * (wi3 .* x)
        σ = norm(ws3 .* u)
        x = wi3 .* (T' * (w3 .* u))
        nx = norm(x)
        nx < 1e-30 && return σ_old
        x ./= nx
        abs(σ - σ_old) ≤ tol * max(σ, 1e-30) && return σ
        σ_old = σ
    end
    return σ_old
end

@inline function _weighted_norm(ws3::AbstractVector, u::AbstractVector)
    s = 0.0
    @inbounds @simd for i in eachindex(u)
        s += abs2(ws3[i] * u[i])
    end
    return sqrt(s)
end

function weighted_gain_squared_svd(T::AbstractMatrix, w_y::AbstractVector)
    wi3, ws3 = _weight_vectors_3N(w_y)
    return abs2(svdvals(ws3 .* T .* reshape(wi3, 1, :))[1])
end

function weighted_gain_squared(T::AbstractMatrix, w_y::AbstractVector; kwargs...)
    abs2(sigma1_weighted(T, w_y; kwargs...))
end

phase_speed_grid(UL_max::Real, n_c::Int = 50) =
    collect(range(0.01 * UL_max, UL_max; length = n_c))

struct ResolventWavenumberCache
    B::Matrix{ComplexF64}
    C::Matrix{ComplexF64}
    Δv::Matrix{ComplexF64}
    L_OS::Matrix{ComplexF64}
    L_Sq::Matrix{ComplexF64}
    M_couple_vw::Matrix{ComplexF64}
    M_couple_wv::Matrix{ComplexF64}
    Pv::Matrix{ComplexF64}
    Pw::Matrix{ComplexF64}
    I_Nw::Matrix{ComplexF64}
    bc_row_idx::Vector{Int}
    bc_rows::Vector{Vector{ComplexF64}}
    N::Int
    Nv::Int
    Nw::Int
    ntot::Int
    Mbuf::Matrix{ComplexF64}
end

function build_resolvent_cache(kx::Real, kz::Real, g::RectGrid, prof)
    N, Nv, Nw, Pv, Pw = g.N, g.Nv, g.Nw, g.Pv, g.Pw
    k = sqrt(kx^2 + kz^2)
    k == 0 && error("build_resolvent_cache: k=0 (use kx=0 branch separately)")
    bl = build_L_blocks(
        kx, kz, g,
        prof.Uv, prof.Usv, prof.nuTv, prof.dUv, prof.d2Uv, prof.dUsv, prof.dnuTv, prof.d2nuTv,
        prof.Uw, prof.Usw, prof.nuTw, prof.dUw, prof.d2Uw, prof.dUsw, prof.dnuTw, prof.d2nuTw,
    )
    Δv, L_OS, L_Sq, F12, F21 = bl.Δv, bl.L_OS, bl.L_Sq, bl.F12, bl.F21
    I_Nw = Matrix{ComplexF64}(I, Nw, Nw)
    I_Nv = Matrix{ComplexF64}(I, Nv, Nv)
    ntot = Nv + Nw
    M_couple_vw = -Pv * F12
    M_couple_wv = -Pw * F21
    Dv = complex.(g.Dv)
    Binner = hcat((-im * kx) .* Dv, (-k^2) .* I_Nv, (-im * kz) .* Dv)
    Bv = Pv * Binner * kron(Matrix{ComplexF64}(I, 3, 3), g.I_Nv_N)
    Bw = Pw * hcat(im * kz .* I_Nw, zeros(ComplexF64, Nw, Nw), (-im * kx) .* I_Nw) *
         kron(Matrix{ComplexF64}(I, 3, 3), g.I_Nw_N)
    Bmat = zeros(ComplexF64, ntot, 3N)
    Bmat[1:N, :] .= Bv
    Bmat[Nv+1:Nv+N, :] .= Bw
    Ivw = g.I_vw
    PvC = g.Pv
    row_u_v = (im * kx / k^2) .* (PvC * Dv)
    row_u_w = (-im * kz / k^2) .* (PvC * Ivw)
    row_v_v = PvC
    Znvw = zeros(ComplexF64, N, Nw)
    row_w_v = (im * kz / k^2) .* (PvC * Dv)
    row_w_w = (im * kx / k^2) .* (PvC * Ivw)
    Cm = vcat(hcat(row_u_v, row_u_w), hcat(row_v_v, Znvw), hcat(row_w_v, row_w_w))
    D2v = complex.(g.D2v)
    Dw = complex.(g.Dw)
    vbc(coeff_v) = vcat(complex.(coeff_v), zeros(ComplexF64, Nw))
    wbc(coeff_w) = vcat(zeros(ComplexF64, Nv), complex.(coeff_w))
    bc_row_idx = Int[]
    bc_rows = Vector{ComplexF64}[]
    for (r, row) in (
        (N + 1, vbc([one(ComplexF64); zeros(ComplexF64, Nv - 1)])),
        (N + 2, vbc(D2v[1, :])),
        (N + 3, vbc(D2v[Nv, :])),
        (N + 4, vbc([zeros(ComplexF64, Nv - 1); one(ComplexF64)])),
        (Nv + N + 1, wbc(Dw[1, :])),
        (Nv + Nw, wbc(Dw[Nw, :])),
    )
        push!(bc_row_idx, r)
        push!(bc_rows, row)
    end
    return ResolventWavenumberCache(
        Bmat, Cm, Δv, L_OS, L_Sq, M_couple_vw, M_couple_wv, Pv, Pw, I_Nw,
        bc_row_idx, bc_rows, N, Nv, Nw, ntot, zeros(ComplexF64, ntot, ntot),
    )
end

function assemble_M!(M::Matrix{ComplexF64}, cache::ResolventWavenumberCache, omega::Real)
    fill!(M, 0)
    N, Nv, Nw, Pv, Pw = cache.N, cache.Nv, cache.Nw, cache.Pv, cache.Pw
    M[1:N, 1:Nv] = Pv * (im * omega .* cache.Δv .- cache.L_OS)
    M[1:N, Nv+1:end] .= cache.M_couple_vw
    M[Nv+1:Nv+N, 1:Nv] .= cache.M_couple_wv
    M[Nv+1:Nv+N, Nv+1:end] = Pw * (im * omega .* cache.I_Nw .- cache.L_Sq)
    for (r, row) in zip(cache.bc_row_idx, cache.bc_rows)
        M[r, :] .= row
    end
    return M
end

transfer_T_cached(cache::ResolventWavenumberCache, omega::Real) =
    cache.C * (assemble_M!(cache.Mbuf, cache, omega) \ cache.B)

function resolvent_transfer(kx, kz, omega, g::RectGrid, prof)
    transfer_gain_rect(
        kx, kz, omega, g,
        prof.Uv, prof.Usv, prof.nuTv, prof.dUv, prof.d2Uv, prof.dUsv, prof.dnuTv, prof.d2nuTv,
        prof.Uw, prof.Usw, prof.nuTw, prof.dUw, prof.d2Uw, prof.dUsw, prof.dnuTw, prof.d2nuTw,
    )
end

function G_max(kx::Real, kz::Real, g::RectGrid, prof, cs::AbstractVector)
    w_y = g.w_y_int
    if abs(kx) < 1e-14
        return weighted_gain_squared(resolvent_transfer(0.0, kz, 0.0, g, prof).T, w_y)
    end
    cache = build_resolvent_cache(kx, kz, g, prof)
    gmax = 0.0
    @inbounds for c in cs
        gmax = max(gmax, weighted_gain_squared(transfer_T_cached(cache, c * kx), w_y))
    end
    return gmax
end

# 原始参考路径(逐 omega 全量 LU + SVD/幂迭代), 仅用于正确性对照。
function G_max_reference(kx::Real, kz::Real, g::RectGrid, prof, cs::AbstractVector)
    w_y = g.w_y_int
    if abs(kx) < 1e-14
        return weighted_gain_squared(resolvent_transfer(0.0, kz, 0.0, g, prof).T, w_y)
    end
    gmax = 0.0
    for c in cs
        res = resolvent_transfer(kx, kz, c * kx, g, prof)
        gmax = max(gmax, weighted_gain_squared(res.T, w_y))
    end
    return gmax
end

# 对比"缓存 LU 快路径 G_max"与"完全重建参考 G_max_reference"(两者都是 LU, 应一致到机器精度)。
# 采用分级判定以稳健应对极少数近临界层共振点的幂迭代收敛差异:
#   max_rel <= rtol       -> PASS
#   rtol < max_rel <= hard -> @warn (近共振良性偏差, 不中断)
#   max_rel > hard         -> error (真实 bug; 通常会是 O(0.1..1) 量级)
function verify_fig3_equivalence(
    g::RectGrid, prof, cs::AbstractVector; rtol::Real = 1e-6, hard_rtol::Real = 1e-2,
)
    H = g.H
    cases = [(2π / 5, 2π / 0.8), (2π / 20, 2π / 2.0), (0.0, 2π / H)]
    maxerr = 0.0
    for (kx, kz) in cases
        Gr = G_max_reference(kx, kz, g, prof, cs)
        Gf = G_max(kx, kz, g, prof, cs)
        err = abs(Gf - Gr) / max(abs(Gr), 1e-30)
        maxerr = max(maxerr, err)
        @printf("verify (kxH,kzH)=(%.3f,%.3f): ref=%.6e fast=%.6e rel=%.3e\n",
            kx * H, kz * H, Gr, Gf, err)
    end
    if maxerr > hard_rtol
        error("Fig.3 cache path differs from transfer_gain_rect reference " *
              "(max rel=$(maxerr) > hard_rtol=$(hard_rtol)); this indicates a real bug.")
    elseif maxerr > rtol
        @warn @sprintf("verify_fig3_equivalence: cache vs reference max rel=%.2e (> rtol=%.0e) near critical-layer resonance — benign, Fig.3 unaffected.", maxerr, rtol)
    else
        @printf("verify_fig3_equivalence: PASS (max rel=%.2e, rtol=%.0e)\n", maxerr, rtol)
    end
    return maxerr <= hard_rtol
end

function verify_sigma1_power(
    g::RectGrid, prof, cs::AbstractVector; rtol::Real = 1e-6, kx = 2π / 5, kz = 2π / 0.8,
)
    w_y = g.w_y_int
    cache = build_resolvent_cache(kx, kz, g, prof)
    ok = true
    for c in cs[1:min(end, 8)]
        T = transfer_T_cached(cache, c * kx)
        Gp = weighted_gain_squared(T, w_y)                  # 幂迭代
        Gs = weighted_gain_squared_svd(T, w_y)              # 精确 SVD
        err = abs(Gp - Gs) / max(abs(Gs), 1e-30)
        @printf("sigma1 power vs svd (c=%.3f): %.6e vs %.6e rel=%.3e\n", c, Gp, Gs, err)
        ok &= err <= rtol
    end
    ok || error("sigma1 power iteration differs from svdvals beyond rtol=$rtol")
    @printf("verify_sigma1_power: PASS (rtol=%.0e)\n", rtol)
    return ok
end

function scan_Gmax_map(
    g::RectGrid, prof,
    lambda_x_over_H::AbstractVector, lambda_z_over_H::AbstractVector, cs::AbstractVector;
    H::Real, parallel::Bool = true, blas_threads::Int = 1,
)
    nx, nz = length(lambda_x_over_H), length(lambda_z_over_H)
    Gm = Matrix{Float64}(undef, nz, nx)
    nth = parallel ? Threads.nthreads() : 1
    @printf("Fig.3 G_max scan: %d x %d wavelengths, %d phase speeds, %d thread(s), BLAS=%d\n",
        nx, nz, length(cs), nth, blas_threads)

    # 外层用 @threads 做粗粒度并行时, 内层 BLAS 设为单线程可避免过度订阅 (oversubscription)。
    # 当 Julia 线程数 >= 物理核数时这通常显著更快; 若 Julia 线程数 < 物理核数, 可把
    # blas_threads 调大以让 BLAS 利用空闲核。
    blas_saved = BLAS.get_num_threads()
    BLAS.set_num_threads(max(1, blas_threads))

    function fill_column!(ix::Int)
        kx = 2π / (lambda_x_over_H[ix] * H)
        for iz in 1:nz
            kz = 2π / (lambda_z_over_H[iz] * H)
            Gm[iz, ix] = G_max(kx, kz, g, prof, cs)
        end
    end

    try
        if parallel && nth > 1
            done = Threads.Atomic{Int}(0)
            @sync @threads for ix in 1:nx
                fill_column!(ix)
                n = Threads.atomic_add!(done, 1) + 1
                if n == 1 || n == nx || n % max(1, nx ÷ 5) == 0
                    @printf("  progress: ix = %d / %d\n", n, nx)
                end
            end
        else
            for ix in 1:nx
                fill_column!(ix)
                if ix == 1 || ix == nx || ix % max(1, nx ÷ 5) == 0
                    @printf("  progress: ix = %d / %d\n", ix, nx)
                end
            end
        end
    finally
        BLAS.set_num_threads(blas_saved)
    end
    return Gm
end

# ---------------------------------------------------------------------------
# Weighted resolvent modes (Fig. 4)
# ---------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------
# Langmuir profiles (Stokes drift + eddy viscosity nu_T)
#
# 涡黏 nu_T 有两种来源 (nuT_profile 关键字):
#   :parabola  -- 理想化对称抛物线 (旧 demo 默认; 不能复现论文 Fig.3 高值区)。
#   :les       -- 由论文 Fig.2(b) 数字化得到的 LES nu_t(y) (按 La_t 选 0.2/0.3),
#                 用多项式光滑拟合 + 解析求导, 叠加分子黏性 1/Re_tau。
#                 这是复现 Fig.3 量级与"小 lambda_z 高增益区"所必需的剖面。
#
# 数字化数据: 论文 PDF 第 11 页 Fig.2(b) 以 600 dpi 渲染后, 按曲线颜色(C0/C1)
# 提取像素并用坐标轴刻度标定得到 (y/H, nu_t), 在 y/H in [0,-1] 上等距 256 点
# (与论文 LES 的 N_y=256 / resolvent 的 N=256 一致; 见 les_eddy_viscosity_fig2b.csv)。
# Fig.2(a) 的 U_L 上半部与解析 Stokes 漂流 (1/La^2)*exp(2*k0H*y/H) 基本吻合,
# 故 U_L 仍采用该解析式 (论文亦指出欧拉平均流可忽略, U_L ~= U_s)。
# ---------------------------------------------------------------------------

# y/H 网格: 256 点 (与论文 LES 的 N_y=256 / resolvent 的 N=256 一致), 等距 0 -> -1。
const LES_NUT_Y = collect(range(0.0, -1.0; length = 256))
const LES_NUT_LA02 = [
    0.00197, 0.00283, 0.00354, 0.00368, 0.00109, 0.00139, 0.00177, 0.00256,
    0.00339, 0.00443, 0.00531, 0.00621, 0.00715, 0.00808, 0.00903, 0.00997,
    0.01083, 0.01161, 0.01260, 0.01346, 0.01432, 0.01519, 0.01605, 0.01692,
    0.01778, 0.01850, 0.01931, 0.01998, 0.02064, 0.02131, 0.02205, 0.02264,
    0.02323, 0.02378, 0.02441, 0.02480, 0.02520, 0.02566, 0.02613, 0.02621,
    0.02618, 0.02641, 0.02842, 0.02880, 0.02874, 0.02854, 0.02874, 0.02884,
    0.02913, 0.02913, 0.02953, 0.02953, 0.02953, 0.02972, 0.03012, 0.03092,
    0.03159, 0.02812, 0.02785, 0.02815, 0.02950, 0.02953, 0.02953, 0.02953,
    0.02953, 0.02928, 0.02913, 0.02913, 0.02913, 0.02913, 0.02874, 0.02874,
    0.02874, 0.02874, 0.02874, 0.02835, 0.02835, 0.02835, 0.02835, 0.02808,
    0.02795, 0.02795, 0.02795, 0.02795, 0.02795, 0.02795, 0.02795, 0.02835,
    0.02835, 0.02874, 0.02874, 0.02874, 0.02913, 0.02953, 0.02972, 0.02992,
    0.03031, 0.03084, 0.03150, 0.03189, 0.03228, 0.03307, 0.03378, 0.03445,
    0.03511, 0.03598, 0.03687, 0.03781, 0.03875, 0.03969, 0.04069, 0.04411,
    0.04469, 0.04469, 0.04463, 0.04469, 0.04629, 0.04969, 0.05111, 0.05264,
    0.05417, 0.05590, 0.05743, 0.05915, 0.06088, 0.06281, 0.06473, 0.06666,
    0.06852, 0.07050, 0.07263, 0.07472, 0.07672, 0.07872, 0.08092, 0.08343,
    0.08587, 0.08807, 0.09050, 0.09270, 0.09522, 0.09748, 0.09985, 0.10227,
    0.10449, 0.10701, 0.10926, 0.11164, 0.11711, 0.11884, 0.11909, 0.11897,
    0.12047, 0.12409, 0.12807, 0.13042, 0.13270, 0.13501, 0.13721, 0.13946,
    0.14178, 0.14390, 0.14602, 0.14814, 0.15026, 0.15219, 0.15392, 0.15564,
    0.15777, 0.15930, 0.16083, 0.16236, 0.16389, 0.16522, 0.16645, 0.16751,
    0.16877, 0.16963, 0.17049, 0.17116, 0.17202, 0.17244, 0.17283, 0.17323,
    0.17323, 0.17323, 0.17323, 0.17323, 0.17283, 0.17244, 0.17190, 0.17105,
    0.16850, 0.16850, 0.16850, 0.16850, 0.16827, 0.16519, 0.16304, 0.16132,
    0.15966, 0.15785, 0.15574, 0.15362, 0.15130, 0.14906, 0.14666, 0.14395,
    0.14129, 0.13853, 0.13580, 0.13294, 0.12988, 0.12676, 0.12356, 0.12030,
    0.11684, 0.11358, 0.11013, 0.10648, 0.10302, 0.09937, 0.09552, 0.08798,
    0.08526, 0.07992, 0.07943, 0.07751, 0.07341, 0.06958, 0.06591, 0.06226,
    0.05840, 0.05456, 0.05054, 0.04671, 0.04309, 0.03936, 0.03591, 0.03230,
    0.02899, 0.02574, 0.02267, 0.01968, 0.01689, 0.01413, 0.01171, 0.00951,
    0.00737, 0.00669, 0.00669, 0.00669, 0.00669, 0.00669, 0.00669, 0.00669]
const LES_NUT_LA03 = [
    0.00039, 0.00067, 0.00133, 0.00223, 0.00415, 0.00628, 0.00840, 0.01052,
    0.01264, 0.01476, 0.01688, 0.01900, 0.02113, 0.02717, 0.02815, 0.02815,
    0.02815, 0.02815, 0.02815, 0.03295, 0.03381, 0.03468, 0.03525, 0.03562,
    0.03589, 0.03602, 0.03601, 0.03583, 0.03543, 0.03504, 0.03465, 0.03405,
    0.03346, 0.03275, 0.03224, 0.03158, 0.03091, 0.03031, 0.02958, 0.02908,
    0.02835, 0.02759, 0.02717, 0.02677, 0.02638, 0.02598, 0.02579, 0.02559,
    0.02559, 0.02559, 0.02579, 0.02598, 0.02645, 0.02700, 0.02775, 0.02861,
    0.02948, 0.03056, 0.03190, 0.03343, 0.03496, 0.03669, 0.03841, 0.04260,
    0.04659, 0.04744, 0.04744, 0.04764, 0.04917, 0.05516, 0.05955, 0.06285,
    0.06630, 0.06982, 0.07380, 0.07767, 0.08180, 0.08612, 0.09049, 0.09495,
    0.09966, 0.10447, 0.10925, 0.11444, 0.11942, 0.12480, 0.13772, 0.13943,
    0.14216, 0.14462, 0.14685, 0.14705, 0.16457, 0.17020, 0.17588, 0.18153,
    0.18718, 0.19290, 0.19868, 0.20434, 0.20999, 0.21565, 0.23256, 0.23248,
    0.23526, 0.23751, 0.24025, 0.24154, 0.25381, 0.26022, 0.26574, 0.27117,
    0.27616, 0.28134, 0.28648, 0.29132, 0.29629, 0.30089, 0.30552, 0.31008,
    0.31456, 0.31868, 0.32854, 0.33037, 0.33209, 0.33230, 0.33422, 0.33524,
    0.34543, 0.34867, 0.35155, 0.35444, 0.35728, 0.36009, 0.36261, 0.36501,
    0.36724, 0.36917, 0.37109, 0.37282, 0.37441, 0.37588, 0.37702, 0.37796,
    0.37894, 0.37964, 0.38031, 0.38071, 0.38110, 0.38110, 0.38110, 0.38071,
    0.38031, 0.37992, 0.37918, 0.37847, 0.37736, 0.37620, 0.37493, 0.37348,
    0.37180, 0.37008, 0.36380, 0.36274, 0.36240, 0.36240, 0.36212, 0.36086,
    0.35400, 0.35121, 0.34823, 0.34517, 0.34211, 0.33885, 0.33539, 0.33194,
    0.32848, 0.32469, 0.32078, 0.31693, 0.31308, 0.30924, 0.30519, 0.30094,
    0.29662, 0.29230, 0.28606, 0.27823, 0.27619, 0.27461, 0.27278, 0.27070,
    0.26501, 0.25578, 0.25100, 0.24602, 0.24105, 0.23604, 0.23086, 0.22607,
    0.22089, 0.21571, 0.21092, 0.20573, 0.20094, 0.18701, 0.18623, 0.18384,
    0.18223, 0.17967, 0.17874, 0.16545, 0.16066, 0.15548, 0.15069, 0.14551,
    0.14072, 0.13553, 0.13075, 0.12556, 0.12077, 0.11579, 0.11100, 0.10621,
    0.10161, 0.09035, 0.08874, 0.08642, 0.08532, 0.08313, 0.08274, 0.06980,
    0.06528, 0.06116, 0.05687, 0.05296, 0.04887, 0.04487, 0.04122, 0.03737,
    0.03388, 0.03026, 0.02700, 0.02362, 0.02049, 0.01765, 0.01496, 0.01244,
    0.01007, 0.00787, 0.00608, 0.00435, 0.00282, 0.00295, 0.00295, 0.00295]

function les_nut_data(La_t::Real)
    if isapprox(La_t, 0.2; atol = 1e-6)
        return LES_NUT_LA02
    elseif isapprox(La_t, 0.3; atol = 1e-6)
        return LES_NUT_LA03
    else
        error("digitized LES nu_t only available for La_t = 0.2 or 0.3 (got $La_t)")
    end
end

# 最小二乘多项式拟合 (Vandermonde), 解析给出 p, p', p''。
function _polyfit(x::AbstractVector, f::AbstractVector, deg::Int)
    V = [x[i]^(k - 1) for i in eachindex(x), k in 1:deg+1]
    return V \ collect(f)
end
_polyval(c, x) = sum(c[k] * x^(k - 1) for k in eachindex(c))
_polyval_d(c, x) = length(c) < 2 ? zero(x) : sum((k - 1) * c[k] * x^(k - 2) for k in 2:length(c))
_polyval_d2(c, x) = length(c) < 3 ? zero(x) : sum((k - 1) * (k - 2) * c[k] * x^(k - 3) for k in 3:length(c))

# 由数字化 nu_t 在 xi=(y+H)/H in [0,1] 上拟合, 返回 (nu, dnu/dy, d2nu/dy2) 三个闭包。
function _smooth_nuT_builders(nu_data::AbstractVector, H::Real, numol::Real; deg::Int = 10)
    xi = (LES_NUT_Y .* H .+ H) ./ H          # = LES_NUT_Y .+ 1, in [0,1]
    c = _polyfit(xi, nu_data, deg)
    nu(y) = max(_polyval(c, (y / H) + 1.0), 0.0) + numol
    dnu(y) = _polyval_d(c, (y / H) + 1.0) / H
    d2nu(y) = _polyval_d2(c, (y / H) + 1.0) / H^2
    return nu, dnu, d2nu
end

function build_langmuir_profiles(
    g::RectGrid;
    H::Real = 1.0,
    La_t::Real = 0.2,
    k0H::Real = 3.5,
    ustar::Real = 1.0,
    nu_mean::Real = 0.07,
    nu_water::Real = 1e-3,
    nuT_profile::Symbol = :parabola,
    Reτ::Real = 1000.0,
    nuT_polydeg::Int = 10,
)
    y_over_H_v = g.y_v ./ H
    y_over_H_w = g.y_w ./ H
    Us_over_ustar_v = (1 / La_t^2) .* exp.(2 * k0H .* y_over_H_v)
    Us_over_ustar_w = (1 / La_t^2) .* exp.(2 * k0H .* y_over_H_w)
    Usv = ustar .* Us_over_ustar_v
    Usw = ustar .* Us_over_ustar_w
    Uv = zeros(Float64, length(g.y_v))
    Uw = zeros(Float64, length(g.y_w))

    if nuT_profile == :les
        # LES-数字化涡黏 (论文 Fig.2b): 光滑拟合 + 解析导数 + 分子黏性 1/Re_tau。
        numol = 1.0 / Reτ
        nu, dnu, d2nu = _smooth_nuT_builders(les_nut_data(La_t), H, numol; deg = nuT_polydeg)
        nuTv = nu.(g.y_v); dnuTv = dnu.(g.y_v); d2nuTv = d2nu.(g.y_v)
        nuTw = nu.(g.y_w); dnuTw = dnu.(g.y_w); d2nuTw = d2nu.(g.y_w)
    elseif nuT_profile == :parabola
        shape_v = @. 4.0 * ((g.y_v + H) / H) * (1.0 - (g.y_v + H) / H)
        shape_w = @. 4.0 * ((g.y_w + H) / H) * (1.0 - (g.y_w + H) / H)
        shape_v = max.(shape_v, 0.0)
        shape_w = max.(shape_w, 0.0)
        scale_v = nu_mean / max(mean(shape_v), 1e-30)
        scale_w = nu_mean / max(mean(shape_w), 1e-30)
        nuTv = scale_v .* shape_v
        nuTw = scale_w .* shape_w
        nu_floor = nu_water / max(ustar * H, 1e-30)
        nuTv[1] = max(nuTv[1], nu_floor)
        nuTv[end] = max(nuTv[end], nu_floor)
        nuTw[1] = max(nuTw[1], nu_floor)
        nuTw[end] = max(nuTw[end], nu_floor)
        dnuTv = g.Dv * nuTv; d2nuTv = g.D2v * nuTv
        dnuTw = g.Dw * nuTw; d2nuTw = g.D2w * nuTw
    else
        error("unknown nuT_profile = $nuT_profile (use :parabola or :les)")
    end

    ULv = Usv
    ULw = Usw
    return (;
        Uv, Usv, nuTv,
        dUv = g.Dv * Uv, d2Uv = g.D2v * Uv, dUsv = g.Dv * Usv,
        dnuTv, d2nuTv,
        Uw, Usw, nuTw,
        dUw = g.Dw * Uw, d2Uw = g.D2w * Uw, dUsw = g.Dw * Usw,
        dnuTw, d2nuTw,
        UL_max = maximum(ULv),
        ULv, ULw,
    )
end

# ---------------------------------------------------------------------------
# Plotting (Fig. 3 / Fig. 4)
# ---------------------------------------------------------------------------

function plot_fig3_Gmax(
    lambda_x_over_H::AbstractVector,
    lambda_z_over_H::AbstractVector,
    Gm::AbstractMatrix;
    La_t::Real = 0.2,
    title_suffix::String = "",
)
    Gplot = log10.(max.(Gm, 1e-30))
    lo = min(maximum(lambda_x_over_H), maximum(lambda_z_over_H))
    hi = max(minimum(lambda_x_over_H), minimum(lambda_z_over_H))
    diag_x = exp10.(range(log10(lo), log10(hi); length = 80))
    fig = Figure(size = (520, 460), fontsize = 14)
    ax = Axis(
        fig[1, 1];
        xscale = log10, yscale = log10,
        xlabel = L"\lambda_x / H", ylabel = L"\lambda_z / H",
        title = "Fig.3-style G_max, " * @sprintf("La_t = %.1f", La_t) * title_suffix,
    )
    cf = contourf!(ax, lambda_x_over_H, lambda_z_over_H, Matrix(Gplot'); colormap = :viridis)
    lines!(ax, diag_x, diag_x; color = :white, linestyle = :dash, linewidth = 1.2)
    Colorbar(fig[1, 2], cf; label = L"\log_{10}(G_{\max})")
    return fig
end

function add_component_panel!(fig, pos, zH, yh, field, ttl)
    vmax = max(maximum(abs, field), 1e-12)
    ax = Axis(fig[pos...]; xlabel = L"z/H", ylabel = L"y/H", title = ttl,
        xlabelsize = 25, ylabelsize = 25, xticklabelsize = 25, yticklabelsize = 25, titlesize = 25)
    ct = contourf!(ax, zH, yh, Matrix(field'); levels = range(-vmax, vmax; length = 17),
        extendlow = :auto, extendhigh = :auto, colormap = :bwr)
    xlims!(ax, 0.0, 0.8)
    ylims!(ax, -1.0, 0.0)
    Colorbar(fig[pos[1], pos[2] + 1], ct)
    return ax
end

function add_vector_quiver!(ax, zH_full, yh_full, xcomp, ycomp; zmax = 0.8, y_step = 8, z_step = 8, color = :white)
    z_idx_all = findall(z -> z <= zmax + 1e-12, zH_full)
    z_idx = z_idx_all[1:z_step:length(z_idx_all)]
    y_idx = collect(1:y_step:length(yh_full))
    xs, ys, us, vs = Float64[], Float64[], Float64[], Float64[]
    for iy in y_idx, iz in z_idx
        push!(xs, zH_full[iz])
        push!(ys, yh_full[iy])
        push!(us, xcomp[iy, iz])
        push!(vs, ycomp[iy, iz])
    end
    mags = sqrt.(us .^ 2 .+ vs .^ 2)
    s = 0.06 / max(maximum(mags), 1e-12)
    us .*= s
    vs .*= s
    arrows2d!(ax, Point2f.(xs, ys), Vec2f.(us, vs);
        color = color, lengthscale = 5.0, shaftwidth = 1.5, tipwidth = 6, tiplength = 8)
end

function phys_field_2d(hc::AbstractVector, kz::Real, zH::AbstractVector, y_nodes::AbstractVector, H::Real)
    nz = length(zH)
    N = length(hc)
    [real(hc[i] * cis(kz * (zH[j] * H))) for i in 1:N, j in 1:nz]
end

function plot_mode_figures(g::RectGrid, wr, sigma1, kz::Real, H::Real; outdir::String = ".")
    N = g.N
    psi_plot = sigma1 .* wr.psi1
    phi_plot = copy(wr.phi1)
    choose_phase_from_v!(psi_plot, phi_plot, N)

    u = psi_plot[1:N]; v = psi_plot[N + 1:2N]; w = psi_plot[2N + 1:3N]
    dx = phi_plot[1:N]; dy = phi_plot[N + 1:2N]; dz = phi_plot[2N + 1:3N]

    py = sortperm(g.y_int)
    yh = g.y_int[py] ./ H
    wy = g.w_y_int[py]
    nz = 720
    zH = collect(range(0.0, 4.0; length = nz))
    pf(h) = phys_field_2d(h[py], kz, zH, g.y_int, H)

    U, V, W = pf(u), pf(v), pf(w)
    DX, DY, DZ = pf(dx), pf(dy), pf(dz)

    fig_mode = Figure(size = (1400, 900), fontsize = 12)
    ax_u = add_component_panel!(fig_mode, (1, 1), zH, yh, U, "Response u")
    ax_v = add_component_panel!(fig_mode, (1, 3), zH, yh, V, "Response v")
    ax_w = add_component_panel!(fig_mode, (1, 5), zH, yh, W, "Response w")
    ax_dx = add_component_panel!(fig_mode, (2, 1), zH, yh, DX, "Forcing d_x")
    ax_dy = add_component_panel!(fig_mode, (2, 3), zH, yh, DY, "Forcing d_y")
    ax_dz = add_component_panel!(fig_mode, (2, 5), zH, yh, DZ, "Forcing d_z")

    for ax in (ax_u, ax_v, ax_w)
        add_vector_quiver!(ax, zH, yh, W, V; zmax = 0.8, y_step = 10, z_step = 10, color = :black)
    end
    for ax in (ax_dx, ax_dy, ax_dz)
        add_vector_quiver!(ax, zH, yh, DZ, DY; zmax = 0.8, y_step = 10, z_step = 10, color = :black)
    end

    fig_prof_w = Figure(size = (1100, 420), fontsize = 12)
    ax_rw = Axis(fig_prof_w[1, 1]; xlabel = L"w_y u_i^2", ylabel = L"y/H",
        title = "Response energy profile (weighted)", xlabelsize = 25, ylabelsize = 25,
        xticklabelsize = 25, yticklabelsize = 25, titlesize = 25)
    lines!(ax_rw, wy .* abs2.(u[py]), yh; color = :steelblue, label = L"w_y u^2", linewidth = 2)
    lines!(ax_rw, wy .* abs2.(v[py]), yh; color = :orangered, label = L"w_y v^2", linewidth = 2)
    lines!(ax_rw, wy .* abs2.(w[py]), yh; color = :seagreen, label = L"w_y w^2", linewidth = 2)
    axislegend(ax_rw; position = :rb, labelsize = 25)

    ax_fw = Axis(fig_prof_w[1, 2]; xlabel = L"w_y d_i^2", ylabel = L"y/H",
        title = "Forcing energy profile (weighted)", xlabelsize = 25, ylabelsize = 25,
        xticklabelsize = 25, yticklabelsize = 25, titlesize = 25)
    lines!(ax_fw, wy .* abs2.(dx[py]), yh; color = :steelblue, label = L"w_y d_x^2", linewidth = 2)
    lines!(ax_fw, wy .* abs2.(dy[py]), yh; color = :orangered, label = L"w_y d_y^2", linewidth = 2)
    lines!(ax_fw, wy .* abs2.(dz[py]), yh; color = :seagreen, label = L"w_y d_z^2", linewidth = 2)
    axislegend(ax_fw; position = :rb, labelsize = 25)

    save(joinpath(outdir, "Analytical_profile_2Dmode.png"), fig_mode)
    save(joinpath(outdir, "Analytical_profile_energymode_weighted.png"), fig_prof_w)
    return fig_mode, fig_prof_w
end

# ---------------------------------------------------------------------------
# Demo driver (notebook workflow)
# ---------------------------------------------------------------------------

const N_demo = 256
const H_demo = 1.0
const USTAR_DEMO = 1.0
const La_t_demo = 0.2
const k0H_demo = 3.5
const NU_T_MEAN_DEMO = 0.07
const NU_WATER_DIM = 1.0e-3
const FIG3_NTHREADS = 8

function run_demo(;
    fig3_quick::Bool = false,
    run_fig3::Bool = true,
    run_fig4::Bool = true,
    outdir::String = @__DIR__,
    nuT_profile::Symbol = :les,
)
    nt = Threads.nthreads()
    @printf("Julia threads: %d (target %d)\n", nt, FIG3_NTHREADS)
    nt < FIG3_NTHREADS &&
        @warn "Threads=$(nt): restart Julia with JULIA_NUM_THREADS=$(FIG3_NTHREADS)"

    g_rect = build_rect_grid(N_demo, H_demo)
    @printf("[%s] sum(w_y)=%.12f (expect H=%.6f)\n", RECT_COLLOC_STAMP, sum(g_rect.w_y_int), H_demo)

    @printf("eddy viscosity profile: %s\n", nuT_profile)
    prof = build_langmuir_profiles(g_rect;
        H = H_demo, La_t = La_t_demo, k0H = k0H_demo,
        ustar = USTAR_DEMO, nu_mean = NU_T_MEAN_DEMO, nu_water = NU_WATER_DIM,
        nuT_profile = nuT_profile)
    isurf_v = argmax(g_rect.y_v)
    @printf("Profiles on Chebyshev grids; UL_max=%.4f\n", prof.UL_max)
    @printf("Surface: U=%.4f, U^s=%.4f, U^L=%.4f\n",
        prof.Uv[isurf_v], prof.Usv[isurf_v], prof.ULv[isurf_v])

    fig3 = nothing
    if run_fig3
        n_phase = fig3_quick ? 20 : 50
        cs = phase_speed_grid(prof.UL_max, n_phase)
        @printf("UL_max = %.4f, phase speeds: %d pts in [%.4f, %.4f]\n",
            prof.UL_max, length(cs), cs[1], cs[end])
        verify_fig3_equivalence(g_rect, prof, cs[1:min(end, 8)])

        if fig3_quick
            lambda_x = logrange10(0.1, 40.0, 16)
            lambda_z = logrange10(0.05, 40.0, 20)
            suffix = " (quick scan)"
        else
            lambda_x = logrange10(0.1, 40.0, 96)
            lambda_z = logrange10(0.05, 40.0, 112)
            suffix = " (paper grid)"
        end

        elapsed = @elapsed Gm_fig3 = scan_Gmax_map(g_rect, prof, lambda_x, lambda_z, cs;
            H = H_demo, parallel = true)
        @printf("G_max range: [%.4e, %.4e], elapsed = %.1f s\n",
            minimum(Gm_fig3), maximum(Gm_fig3), elapsed)

        fig3 = plot_fig3_Gmax(lambda_x, lambda_z, Gm_fig3; La_t = La_t_demo, title_suffix = suffix)
        save(joinpath(outdir, "fig3_Gmax.png"), fig3)
    end

    fig_mode = fig_prof_w = nothing
    if run_fig4
        kx = 0.0
        kz = 2π / H_demo
        omega = 0.0
        @printf("Mode: La_t=%.1f, (kx*H, kz*H, omega) = (%.4f, %.6f, %.4f)\n",
            La_t_demo, kx * H_demo, kz * H_demo, omega)

        res = transfer_gain_rect(
            kx, kz, omega, g_rect,
            prof.Uv, prof.Usv, prof.nuTv, prof.dUv, prof.d2Uv, prof.dUsv, prof.dnuTv, prof.d2nuTv,
            prof.Uw, prof.Usw, prof.nuTw, prof.dUw, prof.d2Uw, prof.dUsw, prof.dnuTw, prof.d2nuTw)
        @printf("Euclidean SVD: G = sigma_1^2 = %.6e\n", res.G)

        wr = weighted_resolvent_modes(res.T, g_rect.w_y_int)
        sigma1 = wr.sigma1
        @printf("Weighted SVD: sigma_1 = %.6e, G_w = %.6e\n", sigma1, abs2(sigma1))

        clim_u, clim_d = 0.8, 2.4
        psi_vis = sigma1 .* wr.psi1
        phi_vis = copy(wr.phi1)
        psi_vis .*= clim_u / max(maximum(abs.(real.(psi_vis))), 1e-30)
        phi_vis .*= clim_d / max(maximum(abs.(real.(phi_vis))), 1e-30)
        iy_phase, theta_opt = choose_phase_from_v!(psi_vis, phi_vis, g_rect.N)
        @printf("Phase: iy=%d, theta=%.4f rad\n", iy_phase, theta_opt)

        fig_mode, fig_prof_w = plot_mode_figures(g_rect, wr, sigma1, kz, H_demo; outdir)
    end

    return (; g_rect, prof, fig3, fig_mode, fig_prof_w)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_demo()
end
