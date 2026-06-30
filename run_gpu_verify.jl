# GPU 批处理后端验证脚本 (在装有 NVIDIA GPU 的机器上运行)。
#
# 运行 (建议开多线程以并行 CPU 端束构造):
#   JULIA_NUM_THREADS=8 julia run_gpu_verify.jl
#   # 或在 REPL 里:  include("run_gpu_verify.jl")
#
# 可选: 用环境变量选卡 (按 nvidia-smi 选空闲的, 例如 H100 通常是 2 号):
#   GPU_ID=2 JULIA_NUM_THREADS=8 julia run_gpu_verify.jl
#
# 依赖同目录下: langmuir_resolvent_cuda.jl 与 langmuir_resolvent_cpu.jl, 以及
# 数字化涡黏 CSV les_eddy_viscosity_fig2b.csv (用 nuT_profile=:les 时)。

include(joinpath(@__DIR__, "langmuir_resolvent_cuda.jl"))
using CUDA, Printf

function main()
    if !CUDA.functional()
        error("CUDA 不可用 (CUDA.functional()==false): 请在装有 NVIDIA 驱动/GPU 的机器上运行。")
    end

    # 选卡
    gpu_id = parse(Int, get(ENV, "GPU_ID", "0"))
    CUDA.device!(gpu_id)
    dev = CUDA.device()
    @printf("使用 GPU #%d: %s\n", gpu_id, CUDA.name(dev))
    @printf("Julia 线程数 = %d (建议 = 物理核数, 用于并行 CPU 端束构造)\n", Threads.nthreads())
    free, total = CUDA.available_memory(), CUDA.total_memory()
    @printf("显存: 可用 %.1f GB / 共 %.1f GB\n", free/2^30, total/2^30)

    # 基流 (与 Fig.3 一致: LES 数字化涡黏)
    H = 1.0
    g = build_rect_grid(256, H)
    prof = build_langmuir_profiles(g; La_t = 0.2, nuT_profile = :les)
    @printf("\nUL_max = %.4f\n", prof.UL_max)

    # 1) 数值正确性: GPU vs CPU
    println("\n================  verify_cuda (GPU vs CPU)  ================")
    cs_v = phase_speed_grid(prof.UL_max, 20)
    verify_cuda(g, prof; cs = cs_v, rtol = 1e-4, maxiter = 300)

    # 2) 计时对比 (小网格)
    println("\n================  bench_cuda (timing)  ================")
    bench_cuda(g, prof; nlx = 8, nlz = 8, ncs = 20, H = H)

    println("\n完成。若 verify_cuda 显示 PASS, 即可用 scan_Gmax_map_cuda 跑完整 Fig.3 扫描:")
    println("    lx = logrange10(0.1,40.0,96); lz = logrange10(0.05,40.0,112)")
    println("    cs = phase_speed_grid(prof.UL_max, 50)")
    println("    Gm = scan_Gmax_map_cuda(g, prof, lx, lz, cs; H=1.0, chunk_pairs=64)")
    return nothing
end

main()
