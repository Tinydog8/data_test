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
using CUDA.CUBLAS: getrf_batched!, getrs_batched!, gemm_strided_batched!

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
    w3::CuVector{Float64}; maxiter::Int = 200,
)
    n3 = size(Tb, 1)
    ncs = size(Tb, 3)
    Xv = CuArray{ComplexF64}(undef, n3, ncs)
    host_init = ComplexF64[cis(2π * i / n3) for i in 1:n3]
    Xv .= CuArray(host_init)                      # 同一初值广播到每列
    Xv ./= sqrt.(sum(abs2, Xv; dims = 1))         # 逐列归一
    Z = similar(Xv); U = similar(Xv); WU = similar(Xv)
    σ = CUDA.zeros(Float64, 1, ncs)
    for _ in 1:maxiter
        @. Z = wi3 * Xv
        gemm_strided_batched!('N', 'N', ComplexF64(1), Tb,
            reshape(Z, n3, 1, ncs), ComplexF64(0), reshape(U, n3, 1, ncs))
        σ .= sqrt.(sum(abs2, ws3 .* U; dims = 1))
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

    # 1) 批量装配: 用 *独立拥有内存* 的 CuMatrix 向量 (不用 view 切片, 兼容老版本
    #    CUDA.jl 的 unsafe_batch 设备指针处理, 避免段错误)。
    Ms = [dA .+ ComplexF64(im * (c * kx)) .* dE for c in cs]   # ncs 个 ntot×ntot

    # 2) 批量 LU: 自分配 pivot 数组并传入 (不依赖返回值顺序)。
    pivots = CuArray{Cint}(undef, ntot, ncs)
    getrf_batched!(Ms, pivots)

    # 3) 批量求解 Xs[j] = M_j^{-1} B  (每个 batch 独立的 B 副本)。
    #    getrs_batched! 签名: (trans, A, B, pivots)。
    Xs = [copy(dB) for _ in 1:ncs]
    getrs_batched!('N', Ms, Xs, pivots)

    # 4) T_j = C * Xs[j]  (写入连续 3D 缓冲, 供批量幂迭代)。
    Tb = CuArray{ComplexF64}(undef, n3, n3, ncs)
    @inbounds for j in 1:ncs
        @views mul!(Tb[:, :, j], dC, Xs[j])
    end

    # 5) 批量加权幂迭代
    σ = _batched_sigma1_gpu(Tb, dwi3, dws3, dw3; maxiter = maxiter)
    return maximum(abs2, σ)
end

# ---------------------------------------------------------------------------
# 微型自检: 在小随机系统上验证 batched LU 求解原语 (getrf_batched!/getrs_batched!)
# 在本机 CUDA.jl 版本上能正常工作。先跑这个; 若它都 segfault/报错, 说明是批处理原语
# 的版本兼容问题, 请把输出贴回。
# ---------------------------------------------------------------------------
function cuda_microtest(; n::Int = 5, ncs::Int = 3, nrhs::Int = 4)
    CUDA.functional() || error("CUDA 不可用。")
    println("cuda_microtest: n=$n, batch=$ncs, nrhs=$nrhs")
    Ahost = [Matrix{ComplexF64}(I, n, n) .+ 0.1 .* randn(ComplexF64, n, n) for _ in 1:ncs]
    Bhost = [randn(ComplexF64, n, nrhs) for _ in 1:ncs]
    Xref = [Ahost[j] \ Bhost[j] for j in 1:ncs]

    Ms = [CuArray(copy(Ahost[j])) for j in 1:ncs]
    Xs = [CuArray(copy(Bhost[j])) for j in 1:ncs]
    pivots = CuArray{Cint}(undef, n, ncs)
    getrf_batched!(Ms, pivots)
    CUDA.synchronize()
    getrs_batched!('N', Ms, Xs, pivots)
    CUDA.synchronize()

    err = maximum(j -> norm(Array(Xs[j]) - Xref[j]) / norm(Xref[j]), 1:ncs)
    @printf("  batched solve max rel err vs CPU = %.3e  -> %s\n", err, err < 1e-8 ? "OK" : "FAIL")
    return err < 1e-8
end

# 设备端权重向量 (与 _weight_vectors_3N 一致)。
function _device_weights(w_y::AbstractVector)
    wi3, ws3 = _weight_vectors_3N(w_y)
    return CuArray(wi3), CuArray(ws3), CuArray(ws3 .* ws3)
end

function G_max_cuda(kx::Real, kz::Real, g::RectGrid, prof, cs::AbstractVector; maxiter::Int = 200)
    w_y = g.w_y_int
    if abs(kx) < 1e-14
        # 流向不变模态: 单点, 直接走 CPU (cheap)
        return weighted_gain_squared(resolvent_transfer(0.0, kz, 0.0, g, prof).T, w_y)
    end
    cache = build_resolvent_cache(kx, kz, g, prof)
    A, E = _build_pencil(cache)
    dwi3, dws3, dw3 = _device_weights(w_y)
    return gmax_cuda_pencil(A, E, cache.B, cache.C, dwi3, dws3, dw3, cs, kx; maxiter = maxiter)
end

# ---------------------------------------------------------------------------
# Fig.3 扫描 (GPU). CPU 端多线程预构造每个 (kx,kz) 的束 (A,E,B,C), 再由 GPU 逐个批处理。
# 分块以控制 CPU 内存; GPU 串行处理 (单 GPU)。
# ---------------------------------------------------------------------------
function scan_Gmax_map_cuda(
    g::RectGrid, prof,
    lambda_x_over_H::AbstractVector, lambda_z_over_H::AbstractVector, cs::AbstractVector;
    H::Real, chunk_pairs::Int = 64, maxiter::Int = 200,
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
                    Gm[iz, ix] = gmax_cuda_pencil(As[s], Es[s], Bs[s], Cs[s],
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
    ncs::Int = 20, H::Real = 1.0)
    cs = phase_speed_grid(prof.UL_max, ncs)
    lx = logrange10(0.3, 30.0, nlx); lz = logrange10(0.1, 20.0, nlz)
    # warmup
    G_max_cuda(2π / 5, 2π / 0.8, g, prof, cs); G_max(2π / 5, 2π / 0.8, g, prof, cs)
    t_gpu = @elapsed Gg = scan_Gmax_map_cuda(g, prof, lx, lz, cs; H = H)
    t_cpu = @elapsed Gc = scan_Gmax_map(g, prof, lx, lz, cs; H = H, blas_threads = 1)
    rel = maximum(abs.(Gg .- Gc) ./ max.(abs.(Gc), 1e-30))
    @printf("bench_cuda %dx%d x %d cs: GPU=%.1fs  CPU=%.1fs  speedup x%.1f  (max rel %.2e)\n",
        nlx, nlz, ncs, t_gpu, t_cpu, t_cpu / t_gpu, rel)
    return (; t_gpu, t_cpu, rel)
end
