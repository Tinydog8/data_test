# GPU (CUDA) 批处理后端 -- Langmuir resolvent Fig.3 G_max 扫描.
#
# 思路 (用户要求的"批处理甜区"):
#   对固定 (kx,kz), 传递矩阵满足线性矩阵束  M(omega) = A + i*omega*E
#   (A = M0, E 仅在两个内部块非零; 已在 CPU 端 build_resolvent_cache 预算好)。
#   于是对该 (kx,kz) 的所有相速度 c (omega = c*kx) 一次性:
#     1. 批量装配   Mb[:,:,j] = A + i*omega_j * E         (GPU 广播)
#     2. 批量 LU    cuBLAS getrf_batched!                  (一次算 ncs 个 LU)
#     3. 批量求解   cuBLAS getrs_batched!  -> Xb = M^{-1}B  (同一 B, ncs 个右端块)
#     4. 形成 T_j = C * Xb[:,:,j]                          (写入 3D 缓冲)
#     5. 批量加权幂迭代 (gemm_strided_batched) -> sigma_1[j]
#   G_max = max_j sigma_1[j]^2。
#
# 适配硬件: A100/A800/H100 等数据中心卡 (FP64 强, 显存大), ComplexF64 精度。
#
# !!! 重要: 本文件需要在装有 NVIDIA GPU 的机器上运行 (CUDA.functional() == true)。
#     开发机无 GPU, 作者未在 GPU 上实测; 请先运行 verify_cuda(...) 与 bench_cuda(...)
#     在你的 GPU 上确认数值正确性与加速比后再用于生产。
#
# 用法:
#   include("langmuir_resolvent_cuda.jl")
#   g    = build_rect_grid(256, 1.0)
#   prof = build_langmuir_profiles(g; La_t=0.2, nuT_profile=:les)
#   verify_cuda(g, prof)                       # GPU vs CPU 正确性
#   Gm = scan_Gmax_map_cuda(g, prof, lx, lz, cs; H=1.0)

include(joinpath(@__DIR__, "langmuir_resolvent_cpu.jl"))

using CUDA
# 两条求解路径:
#   :single      逐 omega 用 CUSOLVER 单系统求解 (CuMatrix 的 \)。最稳, 已在 5.9.6 验证可用。
#   :batched_inv 全程批处理: getrf_strided_batched + getri_strided_batched (批量求逆) +
#                gemm_strided_batched (X=M^{-1}B, T=CX)。完全绕开会段错误的 cublasZgetrsBatched,
#                改用 getrf/getri/gemm 批处理原语。更快, 但需先用 cuda_microtest_batched 确认
#                你的环境下 getrf/getri 批处理不崩 (求逆在近临界层精度略差, 对数坐标 Fig.3 可接受)。
using CUDA.CUBLAS: gemm_strided_batched!, getrf_strided_batched!, getri_strided_batched!

# ---------------------------------------------------------------------------
# 由 CPU 端 cache 构造线性束 (A, E) 以及 B, C (CPU 端, 复数双精度)。
# ---------------------------------------------------------------------------
function _build_pencil(cache::ResolventWavenumberCache)
    ntot, N, Nv, Nw = cache.ntot, cache.N, cache.Nv, cache.Nw
    A = assemble_M_fast!(zeros(ComplexF64, ntot, ntot), cache, 0.0)   # M0
    E = zeros(ComplexF64, ntot, ntot)
    @views E[1:N, 1:Nv] .= cache.OS_om
    @views E[Nv+1:Nv+N, Nv+1:Nv+Nw] .= cache.Sq_om
    return A, E
end

# ---------------------------------------------------------------------------
# 批量加权幂迭代: 对 Tb (n3 x n3 x ncs) 的每个切片求 sigma_1
# of diag(ws3) * T_j * diag(wi3)。固定迭代次数 (GPU 上很便宜, 不做逐列早停)。
# ---------------------------------------------------------------------------
function _batched_sigma1_gpu(
    Tb::CuArray{ComplexF64,3}, wi3::CuVector{Float64}, ws3::CuVector{Float64},
    w3::CuVector{Float64}; maxiter::Int = 100, tol::Real = 1e-7, check_every::Int = 8,
)
    n3 = size(Tb, 1)
    ncs = size(Tb, 3)
    Xv = CuArray{ComplexF64}(undef, n3, ncs)
    host_init = ComplexF64[cis(2π * i / n3) for i in 1:n3]
    Xv .= CuArray(host_init)                      # 同一初值广播到每列
    Xv ./= sqrt.(sum(abs2, Xv; dims = 1))         # 逐列归一
    Z = similar(Xv); U = similar(Xv); WU = similar(Xv)
    σ = CUDA.zeros(Float64, 1, ncs)
    σprev = fill(-1.0, ncs)
    for it in 1:maxiter
        @. Z = wi3 * Xv
        gemm_strided_batched!('N', 'N', ComplexF64(1), Tb,
            reshape(Z, n3, 1, ncs), ComplexF64(0), reshape(U, n3, 1, ncs))
        σ .= sqrt.(sum(abs2, ws3 .* U; dims = 1))
        # 收敛判据 (所有列都稳定才停; 每 check_every 次同步一次以省开销)
        if it % check_every == 0
            σh = vec(Array(σ))
            if maximum(abs.(σh .- σprev) ./ max.(σh, 1e-30)) ≤ tol
                return σh
            end
            σprev = σh
        end
        @. WU = w3 * U
        gemm_strided_batched!('C', 'N', ComplexF64(1), Tb,
            reshape(WU, n3, 1, ncs), ComplexF64(0), reshape(Xv, n3, 1, ncs))
        @. Xv = wi3 * Xv
        Xv ./= sqrt.(sum(abs2, Xv; dims = 1))
    end
    return vec(Array(σ))           # 长度 ncs
end

# ---------------------------------------------------------------------------
# 单个 (kx,kz): 批处理整条 omega 扫描, 返回 G_max。
# 传入预构造的设备常量 (dwi3,dws3,dw3) 以复用。
# ---------------------------------------------------------------------------
function gmax_cuda_pencil(
    A::Matrix{ComplexF64}, E::Matrix{ComplexF64},
    B::Matrix{ComplexF64}, C::Matrix{ComplexF64},
    dwi3, dws3, dw3, cs::AbstractVector, kx::Real; maxiter::Int = 200,
)
    ntot = size(A, 1); n3 = size(B, 2); ncs = length(cs)
    dA = CuArray(A); dE = CuArray(E); dB = CuArray(B); dC = CuArray(C)

    # 1)-3) 逐 omega 单系统求解 (CUSOLVER, 经 CuMatrix 的 \), 形成 T_j 存入连续 3D 缓冲。
    #       绕开有问题的 cublasZgetrsBatched。GPU FP64 极快, ncs 个串行求解仍很便宜。
    Tb = CuArray{ComplexF64}(undef, n3, n3, ncs)
    @inbounds for j in 1:ncs
        Mj = dA .+ ComplexF64(im * (cs[j] * kx)) .* dE          # ntot×ntot
        Xj = Mj \ dB                                            # CUSOLVER 单系统 LU 解, ntot×n3
        @views mul!(Tb[:, :, j], dC, Xj)                        # T_j = C X_j
    end

    # 4) 批量加权幂迭代 (gemm_strided_batched, 与 getrsBatched 无关)
    σ = _batched_sigma1_gpu(Tb, dwi3, dws3, dw3; maxiter = maxiter)
    return maximum(abs2, σ)
end

# 全程批处理版: getrf + getri (批量求逆) + gemm 批处理, 绕开 cublasZgetrsBatched。
function gmax_cuda_pencil_batched(
    A::Matrix{ComplexF64}, E::Matrix{ComplexF64},
    B::Matrix{ComplexF64}, C::Matrix{ComplexF64},
    dwi3, dws3, dw3, cs::AbstractVector, kx::Real; maxiter::Int = 100,
)
    ntot = size(A, 1); n3 = size(B, 2); ncs = length(cs)
    dA = CuArray(A); dE = CuArray(E); dB = CuArray(B); dC = CuArray(C)
    iω = CuArray(ComplexF64[im * (c * kx) for c in cs])

    # 1) 批量装配 + 批量 LU
    Mb = dA .+ reshape(iω, 1, 1, ncs) .* dE                     # (ntot,ntot,ncs)
    pivots = CuArray{Cint}(undef, ntot, ncs)
    getrf_strided_batched!(Mb, pivots)

    # 2) 批量求逆 Minv = M^{-1}
    Minv = CuArray{ComplexF64}(undef, ntot, ntot, ncs)
    getri_strided_batched!(Mb, Minv, pivots)

    # 3) X = Minv * B (批量 gemm, B 广播到每个 batch)
    Bb = CuArray{ComplexF64}(undef, ntot, n3, ncs); Bb .= reshape(dB, ntot, n3, 1)
    Xb = CuArray{ComplexF64}(undef, ntot, n3, ncs)
    gemm_strided_batched!('N', 'N', ComplexF64(1), Minv, Bb, ComplexF64(0), Xb)

    # 4) T = C * X (批量 gemm, C 广播到每个 batch)
    Cb = CuArray{ComplexF64}(undef, n3, ntot, ncs); Cb .= reshape(dC, n3, ntot, 1)
    Tb = CuArray{ComplexF64}(undef, n3, n3, ncs)
    gemm_strided_batched!('N', 'N', ComplexF64(1), Cb, Xb, ComplexF64(0), Tb)

    # 5) 批量加权幂迭代
    σ = _batched_sigma1_gpu(Tb, dwi3, dws3, dw3; maxiter = maxiter)
    return maximum(abs2, σ)
end

# ---------------------------------------------------------------------------
# 微型自检: 在小随机系统上验证 batched LU 求解原语 (getrf_batched!/getrs_batched!)
# 在本机 CUDA.jl 版本上能正常工作。先跑这个; 若它都 segfault/报错, 说明是批处理原语
# 的版本兼容问题, 请把输出贴回。
# ---------------------------------------------------------------------------
function cuda_microtest(; n::Int = 256, ncs::Int = 4, nrhs::Int = 200)
    CUDA.functional() || error("CUDA 不可用。")
    println("cuda_microtest: n=$n, batch=$ncs, nrhs=$nrhs")

    # 原语 1: 单系统 CUSOLVER 求解 (\) —— gmax_cuda_pencil 实际使用的求解
    A1 = Matrix{ComplexF64}(I, n, n) .+ 0.1 .* randn(ComplexF64, n, n)
    B1 = randn(ComplexF64, n, nrhs)
    Xref = A1 \ B1
    Xg = CuArray(A1) \ CuArray(B1)
    CUDA.synchronize()
    err1 = norm(Array(Xg) - Xref) / norm(Xref)
    @printf("  [1] single-system CuMatrix \\ : rel err = %.3e  -> %s\n", err1, err1 < 1e-8 ? "OK" : "FAIL")

    # 原语 2: 批量 gemm (gemm_strided_batched!) —— 幂迭代实际使用
    m = 64
    P = randn(ComplexF64, m, m, ncs); Q = randn(ComplexF64, m, 1, ncs)
    dP = CuArray(P); dQ = CuArray(Q); dR = CuArray{ComplexF64}(undef, m, 1, ncs)
    gemm_strided_batched!('N', 'N', ComplexF64(1), dP, dQ, ComplexF64(0), dR)
    CUDA.synchronize()
    Rh = Array(dR)
    err2 = maximum(j -> norm(Rh[:, :, j] - P[:, :, j] * Q[:, :, j]) / norm(P[:, :, j] * Q[:, :, j]), 1:ncs)
    @printf("  [2] gemm_strided_batched!    : rel err = %.3e  -> %s\n", err2, err2 < 1e-8 ? "OK" : "FAIL")

    ok = err1 < 1e-8 && err2 < 1e-8
    println("  microtest ", ok ? "PASS" : "FAIL")
    return ok
end

# 设备端权重向量 (与 _weight_vectors_3N 一致)。
function _device_weights(w_y::AbstractVector)
    wi3, ws3 = _weight_vectors_3N(w_y)
    return CuArray(wi3), CuArray(ws3), CuArray(ws3 .* ws3)
end

_pencil_gain(solve::Symbol, args...; kwargs...) =
    solve === :batched_inv ? gmax_cuda_pencil_batched(args...; kwargs...) :
    gmax_cuda_pencil(args...; kwargs...)

function G_max_cuda(kx::Real, kz::Real, g::RectGrid, prof, cs::AbstractVector;
    maxiter::Int = 100, solve::Symbol = :single)
    w_y = g.w_y_int
    if abs(kx) < 1e-14
        # 流向不变模态: 单点, 直接走 CPU (cheap)
        return weighted_gain_squared(resolvent_transfer(0.0, kz, 0.0, g, prof).T, w_y)
    end
    cache = build_resolvent_cache(kx, kz, g, prof)
    A, E = _build_pencil(cache)
    dwi3, dws3, dw3 = _device_weights(w_y)
    return _pencil_gain(solve, A, E, cache.B, cache.C, dwi3, dws3, dw3, cs, kx; maxiter = maxiter)
end

# ---------------------------------------------------------------------------
# Fig.3 扫描 (GPU). CPU 端多线程预构造每个 (kx,kz) 的束 (A,E,B,C), 再由 GPU 逐个批处理。
# 分块以控制 CPU 内存; GPU 串行处理 (单 GPU)。
# ---------------------------------------------------------------------------
function scan_Gmax_map_cuda(
    g::RectGrid, prof,
    lambda_x_over_H::AbstractVector, lambda_z_over_H::AbstractVector, cs::AbstractVector;
    H::Real, chunk_pairs::Int = 64, maxiter::Int = 100, solve::Symbol = :single,
)
    CUDA.functional() || error("CUDA 不可用 (CUDA.functional()==false); 需在装有 NVIDIA GPU 的机器上运行。")
    nx, nz = length(lambda_x_over_H), length(lambda_z_over_H)
    Gm = Matrix{Float64}(undef, nz, nx)
    dwi3, dws3, dw3 = _device_weights(g.w_y_int)
    @printf("Fig.3 G_max scan (CUDA): %d x %d wavelengths, %d phase speeds, GPU=%s\n",
        nx, nz, length(cs), CUDA.name(CUDA.device()))

    # 把 (kx,kz) 线性化为 pair 列表
    pairs = [(ix, iz) for ix in 1:nx for iz in 1:nz]
    npairs = length(pairs)
    blas_saved = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    try
        idx = 1
        while idx <= npairs
            hi = min(idx + chunk_pairs - 1, npairs)
            sub = pairs[idx:hi]
            ns = length(sub)
            # CPU 多线程构造束 (A,E,B,C) 与 kx 列表
            As = Vector{Matrix{ComplexF64}}(undef, ns)
            Es = Vector{Matrix{ComplexF64}}(undef, ns)
            Bs = Vector{Matrix{ComplexF64}}(undef, ns)
            Cs = Vector{Matrix{ComplexF64}}(undef, ns)
            kxs = Vector{Float64}(undef, ns)
            kx0_flag = falses(ns)
            @threads for s in 1:ns
                ix, iz = sub[s]
                kx = 2π / (lambda_x_over_H[ix] * H)
                kz = 2π / (lambda_z_over_H[iz] * H)
                kxs[s] = kx
                if abs(kx) < 1e-14
                    kx0_flag[s] = true
                else
                    cache = build_resolvent_cache(kx, kz, g, prof)
                    A, E = _build_pencil(cache)
                    As[s] = A; Es[s] = E; Bs[s] = cache.B; Cs[s] = cache.C
                end
            end
            # GPU 串行处理本 chunk
            for s in 1:ns
                ix, iz = sub[s]
                if kx0_flag[s]
                    kz = 2π / (lambda_z_over_H[iz] * H)
                    Gm[iz, ix] = weighted_gain_squared(resolvent_transfer(0.0, kz, 0.0, g, prof).T, g.w_y_int)
                else
                    Gm[iz, ix] = _pencil_gain(solve, As[s], Es[s], Bs[s], Cs[s],
                        dwi3, dws3, dw3, cs, kxs[s]; maxiter = maxiter)
                end
            end
            @printf("  progress: %d / %d pairs\n", hi, npairs)
            idx = hi + 1
        end
    finally
        BLAS.set_num_threads(blas_saved)
    end
    return Gm
end

# ---------------------------------------------------------------------------
# 校验: GPU vs CPU 在若干 (kx,kz) 上对比 G_max。
# ---------------------------------------------------------------------------
function verify_cuda(g::RectGrid, prof; cs::AbstractVector = phase_speed_grid(prof.UL_max, 20),
    rtol::Real = 1e-4, maxiter::Int = 300)
    CUDA.functional() || error("CUDA 不可用; 需在 GPU 机器上运行。")
    H = g.H
    cases = [(2π / 5, 2π / 0.8), (2π / 20, 2π / 2.0), (2π / 3, 2π / 0.5)]
    ok = true
    for (kx, kz) in cases
        Gc = G_max(kx, kz, g, prof, cs)
        Gg = G_max_cuda(kx, kz, g, prof, cs; maxiter = maxiter)
        err = abs(Gg - Gc) / max(abs(Gc), 1e-30)
        @printf("verify_cuda (kxH,kzH)=(%.3f,%.3f): cpu=%.6e gpu=%.6e rel=%.3e\n",
            kx * H, kz * H, Gc, Gg, err)
        ok &= err <= rtol
    end
    ok || error("CUDA G_max differs from CPU beyond rtol=$rtol")
    @printf("verify_cuda: PASS (rtol=%.0e)\n", rtol)
    return ok
end

# ---------------------------------------------------------------------------
# 计时对比: 一小块网格上 CPU vs GPU 扫描。
# ---------------------------------------------------------------------------
function bench_cuda(g::RectGrid, prof; nlx::Int = 8, nlz::Int = 8,
    ncs::Int = 50, H::Real = 1.0, solve::Symbol = :single)
    cs = phase_speed_grid(prof.UL_max, ncs)
    lx = logrange10(0.3, 30.0, nlx); lz = logrange10(0.1, 20.0, nlz)
    npairs = nlx * nlz
    # warmup
    G_max_cuda(2π / 5, 2π / 0.8, g, prof, cs; solve = solve); G_max(2π / 5, 2π / 0.8, g, prof, cs)

    # 诊断: 仅 CPU 端束构造 (并行) 的耗时, 衡量 GPU 版的"地板"
    blas_saved = BLAS.get_num_threads(); BLAS.set_num_threads(1)
    t_build = @elapsed begin
        @threads for ix in 1:nlx
            for iz in 1:nlz
                kx = 2π / (lx[ix] * H); kz = 2π / (lz[iz] * H)
                cache = build_resolvent_cache(kx, kz, g, prof)
                _build_pencil(cache)
            end
        end
    end
    BLAS.set_num_threads(blas_saved)

    t_gpu = @elapsed Gg = scan_Gmax_map_cuda(g, prof, lx, lz, cs; H = H, solve = solve)
    t_cpu = @elapsed Gc = scan_Gmax_map(g, prof, lx, lz, cs; H = H, blas_threads = 1)
    rel = maximum(abs.(Gg .- Gc) ./ max.(abs.(Gc), 1e-30))
    @printf("\nbench_cuda %dx%d=%d pairs x %d cs  (N=%d, %d threads, solve=%s):\n",
        nlx, nlz, npairs, ncs, g.N, Threads.nthreads(), solve)
    @printf("  CPU 束构造(并行, GPU 版地板) : %.1f s\n", t_build)
    @printf("  GPU 全程                     : %.1f s\n", t_gpu)
    @printf("  CPU 全程                     : %.1f s\n", t_cpu)
    @printf("  speedup x%.2f   (max rel %.2e)\n", t_cpu / t_gpu, rel)
    return (; t_build, t_gpu, t_cpu, rel)
end

# ---------------------------------------------------------------------------
# 微型自检 (批处理路径): 验证 getrf_strided_batched + getri_strided_batched + gemm 在本机
# 不崩、且数值正确。若 PASS, 可用 solve=:batched_inv 获得全程批处理加速。
# ---------------------------------------------------------------------------
function cuda_microtest_batched(; n::Int = 256, ncs::Int = 8)
    CUDA.functional() || error("CUDA 不可用。")
    println("cuda_microtest_batched: n=$n, batch=$ncs (getrf+getri+gemm strided batched)")
    Ah = [Matrix{ComplexF64}(I, n, n) .+ 0.1 .* randn(ComplexF64, n, n) for _ in 1:ncs]
    Mb = CuArray{ComplexF64}(undef, n, n, ncs)
    for j in 1:ncs
        @views Mb[:, :, j] .= CuArray(Ah[j])
    end
    pivots = CuArray{Cint}(undef, n, ncs)
    getrf_strided_batched!(Mb, pivots); CUDA.synchronize()
    Minv = CuArray{ComplexF64}(undef, n, n, ncs)
    getri_strided_batched!(Mb, Minv, pivots); CUDA.synchronize()
    Mih = Array(Minv)
    err = maximum(j -> norm(Mih[:, :, j] - inv(Ah[j])) / norm(inv(Ah[j])), 1:ncs)
    ok = err < 1e-7
    @printf("  getrf+getri batched inverse: rel err = %.3e -> %s\n", err, ok ? "OK" : "FAIL")
    println("  -> 若 OK 且不段错误, 可用 solve=:batched_inv 加速。")
    return ok
end
