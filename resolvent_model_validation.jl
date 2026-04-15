using LinearAlgebra
using Printf
using DelimitedFiles

required_symbols = [:build_rect_grid, :build_L_blocks, :build_M_B_C_rect, :transfer_gain_rect]
missing_symbols = [s for s in required_symbols if !isdefined(Main, s)]
isempty(missing_symbols) || error("请先运行前面的模型定义；缺少符号: $(join(string.(missing_symbols), \", \"))")

relerr(a, b) = norm(a .- b) / max(norm(b), 1e-30)

function report_test(name::AbstractString, value::Real, tol::Real)
    ok = value <= tol
    @printf("%-40s : %-4s value = %.3e   tol = %.1e\n", name, ok ? "PASS" : "FAIL", value, tol)
    return ok
end

function constant_case_data(g; ν0 = 0.05)
    Uv = zeros(Float64, length(g.y_v))
    Usv = zeros(Float64, length(g.y_v))
    nuTv = fill(Float64(ν0), length(g.y_v))
    dUv = zeros(Float64, length(g.y_v))
    d2Uv = zeros(Float64, length(g.y_v))
    dUsv = zeros(Float64, length(g.y_v))
    dnuTv = zeros(Float64, length(g.y_v))
    d2nuTv = zeros(Float64, length(g.y_v))

    Uw = zeros(Float64, length(g.y_w))
    Usw = zeros(Float64, length(g.y_w))
    nuTw = fill(Float64(ν0), length(g.y_w))
    dUw = zeros(Float64, length(g.y_w))
    d2Uw = zeros(Float64, length(g.y_w))
    dUsw = zeros(Float64, length(g.y_w))
    dnuTw = zeros(Float64, length(g.y_w))
    d2nuTw = zeros(Float64, length(g.y_w))

    return (Uv, Usv, nuTv, dUv, d2Uv, dUsv, dnuTv, d2nuTv,
            Uw, Usw, nuTw, dUw, d2Uw, dUsw, dnuTw, d2nuTw)
end

function linear_shear_case_data(g; S = 1.25, ν0 = 0.05)
    Uv = S .* g.y_v
    Usv = zeros(Float64, length(g.y_v))
    nuTv = fill(Float64(ν0), length(g.y_v))
    dUv = fill(Float64(S), length(g.y_v))
    d2Uv = zeros(Float64, length(g.y_v))
    dUsv = zeros(Float64, length(g.y_v))
    dnuTv = zeros(Float64, length(g.y_v))
    d2nuTv = zeros(Float64, length(g.y_v))

    Uw = S .* g.y_w
    Usw = zeros(Float64, length(g.y_w))
    nuTw = fill(Float64(ν0), length(g.y_w))
    dUw = fill(Float64(S), length(g.y_w))
    d2Uw = zeros(Float64, length(g.y_w))
    dUsw = zeros(Float64, length(g.y_w))
    dnuTw = zeros(Float64, length(g.y_w))
    d2nuTw = zeros(Float64, length(g.y_w))

    return (Uv, Usv, nuTv, dUv, d2Uv, dUsv, dnuTv, d2nuTv,
            Uw, Usw, nuTw, dUw, d2Uw, dUsw, dnuTw, d2nuTw)
end

function stokes_shear_case_data(g; S = 1.25, ν0 = 0.05)
    Uv = zeros(Float64, length(g.y_v))
    Usv = S .* g.y_v
    nuTv = fill(Float64(ν0), length(g.y_v))
    dUv = zeros(Float64, length(g.y_v))
    d2Uv = zeros(Float64, length(g.y_v))
    dUsv = fill(Float64(S), length(g.y_v))
    dnuTv = zeros(Float64, length(g.y_v))
    d2nuTv = zeros(Float64, length(g.y_v))

    Uw = zeros(Float64, length(g.y_w))
    Usw = S .* g.y_w
    nuTw = fill(Float64(ν0), length(g.y_w))
    dUw = zeros(Float64, length(g.y_w))
    d2Uw = zeros(Float64, length(g.y_w))
    dUsw = fill(Float64(S), length(g.y_w))
    dnuTw = zeros(Float64, length(g.y_w))
    d2nuTw = zeros(Float64, length(g.y_w))

    return (Uv, Usv, nuTv, dUv, d2Uv, dUsv, dnuTv, d2nuTv,
            Uw, Usw, nuTw, dUw, d2Uw, dUsw, dnuTw, d2nuTw)
end

function run_internal_validation(; H = 1.0, N = 128, kx = 1.3, kz = 2.1, ω = 0.7, ν0 = 0.05, S = 1.25, nmode = 2)
    g = build_rect_grid(N, H)
    k2 = kx^2 + kz^2
    Dv = complex.(g.Dv)
    pass = Bool[]

    @printf("内部验证开始: H=%.3f, N=%d, kx=%.3f, kz=%.3f, ω=%.3f, ν0=%.3f, S=%.3f\n\n",
        H, N, kx, kz, ω, ν0, S)

    # Case A: U = 0, U^s = 0, ν_T = 常数；耦合块应为 0
    caseA = constant_case_data(g; ν0 = ν0)
    blA = build_L_blocks(kx, kz, g, caseA...)
    push!(pass, report_test("Case A: norm(F12)", norm(blA.F12), 1e-12))
    push!(pass, report_test("Case A: norm(F21)", norm(blA.F21), 1e-12))

    # 常系数情形下，用满足边界条件的解析模态验证 L_OS 与 L_Sq
    s_v = (g.y_v .+ H) ./ H
    s_w = (g.y_w .+ H) ./ H
    vn = complex.(sin.(nmode * pi .* s_v))
    wn = complex.(cos.(nmode * pi .* s_w))
    α2 = k2 + (nmode * pi / H)^2
    λOS = ν0 * α2^2
    λSq = -ν0 * α2
    push!(pass, report_test("Case A: L_OS 本征残差", relerr(blA.L_OS * vn, λOS .* vn), 1e-8))
    push!(pass, report_test("Case A: L_Sq 本征残差", relerr(blA.L_Sq * wn, λSq .* wn), 1e-8))

    # B / C 结构验证
    M, Bmat, Cmat = build_M_B_C_rect(kx, kz, ω, g, caseA...)
    s_int = (g.y_int .+ H) ./ H
    dx = complex.(cos.(pi .* s_int))
    dy = complex.(sin.(2pi .* s_int))
    dz = complex.(cos.(3pi .* s_int))
    d = vcat(dx, dy, dz)
    q = Bmat * d

    dx_v = g.I_Nv_N * dx
    dy_v = g.I_Nv_N * dy
    dz_v = g.I_Nv_N * dz
    dx_w = g.I_Nw_N * dx
    dz_w = g.I_Nw_N * dz

    qv_full = (-im * kx) .* (Dv * dx_v) .- k2 .* dy_v .- (im * kz) .* (Dv * dz_v)
    qω_full = (im * kz) .* dx_w .- (im * kx) .* dz_w
    qv_expected = g.Pv * qv_full
    qω_expected = g.Pw * qω_full
    push!(pass, report_test("B: v 方程 interior 块", relerr(q[1:g.N], qv_expected), 1e-10))
    push!(pass, report_test("B: ω 方程 interior 块", relerr(q[g.Nv+1:g.Nv+g.N], qω_expected), 1e-10))

    ξ = vcat(vn, wn)
    uhat = Cmat * ξ
    u_expected = (im * kx / k2) .* (g.Pv * (Dv * vn)) .- (im * kz / k2) .* (g.Pv * (g.I_vw * wn))
    v_expected = g.Pv * vn
    w_expected = (im * kz / k2) .* (g.Pv * (Dv * vn)) .+ (im * kx / k2) .* (g.Pv * (g.I_vw * wn))
    push!(pass, report_test("C: u 重构", relerr(uhat[1:g.N], u_expected), 1e-10))
    push!(pass, report_test("C: v 重构", relerr(uhat[g.N+1:2g.N], v_expected), 1e-10))
    push!(pass, report_test("C: w 重构", relerr(uhat[2g.N+1:3g.N], w_expected), 1e-10))

    continuity_resid = norm((im * kx) .* uhat[1:g.N] .+ (g.Pv * (Dv * vn)) .+ (im * kz) .* uhat[2g.N+1:3g.N]) / max(norm(uhat), 1e-30)
    vort_expected = g.Pv * (g.I_vw * wn)
    vort_resid = relerr((im * kz) .* uhat[1:g.N] .- (im * kx) .* uhat[2g.N+1:3g.N], vort_expected)
    push!(pass, report_test("C: 连续性残差", continuity_resid, 1e-10))
    push!(pass, report_test("C: ω_y 定义残差", vort_resid, 1e-10))

    rhs = Bmat * d
    x = M \ rhs
    solve_resid = norm(M * x - rhs) / max(norm(rhs), 1e-30)
    push!(pass, report_test("线性求解残差", solve_resid, 1e-12))

    # Case B: U = S y, U^s = 0，应只保留 F21 = -ikz U'
    caseB = linear_shear_case_data(g; S = S, ν0 = ν0)
    blB = build_L_blocks(kx, kz, g, caseB...)
    F21_expected = (-im * kz * S) .* g.I_wv
    push!(pass, report_test("Case B: F12 仍为 0", norm(blB.F12), 1e-12))
    push!(pass, report_test("Case B: F21 与 -ikz*U' 一致", relerr(blB.F21, F21_expected), 1e-12))

    # Case C: U = 0, U^s = S y，应只保留 F12 = -ikz U_s'
    caseC = stokes_shear_case_data(g; S = S, ν0 = ν0)
    blC = build_L_blocks(kx, kz, g, caseC...)
    F12_expected = (-im * kz * S) .* g.I_vw
    push!(pass, report_test("Case C: F21 仍为 0", norm(blC.F21), 1e-12))
    push!(pass, report_test("Case C: F12 与 -ikz*U_s' 一致", relerr(blC.F12, F12_expected), 1e-12))

    # 简单网格收敛
    sigma_list = Float64[]
    for Ntest in (64, 128, 256)
        gt = build_rect_grid(Ntest, H)
        caset = constant_case_data(gt; ν0 = ν0)
        rest = transfer_gain_rect(kx, kz, ω, gt, caset...)
        push!(sigma_list, Float64(real(rest.s1)))
    end
    conv23 = abs(sigma_list[3] - sigma_list[2]) / max(abs(sigma_list[3]), 1e-30)
    @printf("\nσ₁(N=64,128,256) = [%.6e, %.6e, %.6e]\n", sigma_list...)
    push!(pass, report_test("网格收敛: |σ₁(256)-σ₁(128)|/|σ₁(256)|", conv23, 5e-3))

    n_pass = count(identity, pass)
    @printf("\n验证完成: %d / %d 项通过。\n", n_pass, length(pass))
    return (; all_pass = all(pass), pass, sigma_list)
end

function run_analytic_u_pattern_demo(;
    H = 1.0,
    N = 128,
    kz = 2pi,
    ν0 = 0.05,
    nmode = 1,
    nz = 128,
    z_periods = 2.0,
    output_csv::Union{Nothing, AbstractString} = nothing,
)
    kx = 0.0
    ω = 0.0
    g = build_rect_grid(N, H)
    caseA = constant_case_data(g; ν0 = ν0)
    M, Bmat, Cmat = build_M_B_C_rect(kx, kz, ω, g, caseA...)

    s_int = (g.y_int .+ H) ./ H
    dx = complex.(cos.(nmode * pi .* s_int))
    dy = zeros(ComplexF64, g.N)
    dz = zeros(ComplexF64, g.N)
    d = vcat(dx, dy, dz)

    rhs = Bmat * d
    ξ = M \ rhs
    uhat = Cmat * ξ
    u_num_amp = real.(uhat[1:g.N])

    α2 = kz^2 + (nmode * pi / H)^2
    u_exact_amp = cos.(nmode * pi .* s_int) ./ (ν0 * α2)
    amp_relerr = relerr(u_num_amp, u_exact_amp)
    imag_leak = maximum(abs.(imag.(uhat[1:g.N])))

    zH = collect(range(0.0, z_periods; length = nz))
    u_num_pattern = [u_num_amp[i] * cos(kz * H * zH[j]) for i in 1:g.N, j in 1:nz]
    u_exact_pattern = [u_exact_amp[i] * cos(kz * H * zH[j]) for i in 1:g.N, j in 1:nz]
    pattern_relerr = relerr(vec(u_num_pattern), vec(u_exact_pattern))

    @printf("\n典型算例: 常系数 + kx=0 + 仅 d_x 强迫\n")
    @printf("参数: H=%.3f, N=%d, kz=%.6f, ν0=%.4f, nmode=%d\n", H, N, kz, ν0, nmode)
    @printf("解析预期: u(y,z) = cos(%dπ(y+H)/H) * cos(kz z) / (ν0 * α²), α² = kz² + (%dπ/H)²\n", nmode, nmode)
    @printf("u 幅值相对误差      = %.6e\n", amp_relerr)
    @printf("u pattern 相对误差 = %.6e\n", pattern_relerr)
    @printf("u 数值解虚部泄漏最大值 = %.6e\n", imag_leak)

    sample_ids = unique(round.(Int, range(1, g.N; length = 5)))
    println("\n样本点对比: y/H, u_num_amp, u_exact_amp")
    for i in sample_ids
        @printf("% .6f   % .6e   % .6e\n", g.y_int[i] / H, u_num_amp[i], u_exact_amp[i])
    end

    if output_csv !== nothing
        rows = Matrix{Float64}(undef, g.N * nz, 4)
        row = 1
        for i in 1:g.N, j in 1:nz
            rows[row, 1] = g.y_int[i] / H
            rows[row, 2] = zH[j]
            rows[row, 3] = u_num_pattern[i, j]
            rows[row, 4] = u_exact_pattern[i, j]
            row += 1
        end
        writedlm(output_csv, rows, ',')
        @printf("\n已写出 pattern CSV: %s\n", output_csv)
        println("列顺序: y_over_H, z_over_H, u_num, u_exact")
    end

    return (;
        y_over_H = g.y_int ./ H,
        z_over_H = zH,
        u_num_amp,
        u_exact_amp,
        u_num_pattern,
        u_exact_pattern,
        amp_relerr,
        pattern_relerr,
        imag_leak,
    )
end

function save_analytic_u_pattern_plot(
    demo_result;
    output_png::AbstractString = "analytic_u_pattern_demo.png",
    title_prefix::AbstractString = "Analytic u-pattern demo",
)
    if !isnothing(Base.find_package("CairoMakie"))
        @eval using CairoMakie
    else
        error("未找到 CairoMakie；如需直接绘图，请先安装或在已有绘图库环境中运行。")
    end

    z = demo_result.z_over_H
    y = demo_result.y_over_H
    unum = demo_result.u_num_pattern
    uexact = demo_result.u_exact_pattern
    uerr = unum .- uexact

    umax = max(maximum(abs.(unum)), maximum(abs.(uexact)), 1e-12)
    emax = max(maximum(abs.(uerr)), 1e-12)

    fig = Figure(size = (1200, 420), fontsize = 13)

    ax1 = Axis(fig[1, 1], xlabel = "z/H", ylabel = "y/H", title = title_prefix * " - numerical")
    hm1 = heatmap!(ax1, z, y, transpose(unum); colormap = :balance, colorrange = (-umax, umax))
    Colorbar(fig[1, 2], hm1; label = "u_num")

    ax2 = Axis(fig[1, 3], xlabel = "z/H", ylabel = "y/H", title = title_prefix * " - analytic")
    hm2 = heatmap!(ax2, z, y, transpose(uexact); colormap = :balance, colorrange = (-umax, umax))
    Colorbar(fig[1, 4], hm2; label = "u_exact")

    ax3 = Axis(fig[1, 5], xlabel = "z/H", ylabel = "y/H", title = title_prefix * " - error")
    hm3 = heatmap!(ax3, z, y, transpose(uerr); colormap = :balance, colorrange = (-emax, emax))
    Colorbar(fig[1, 6], hm3; label = "u_num - u_exact")

    save(output_png, fig)
    @printf("已保存二维 pattern 图: %s\n", output_png)
    return fig
end

validation_result = run_internal_validation()
println("\nvalidation_result = ", validation_result)

demo_result = run_analytic_u_pattern_demo()
println("\ndemo_result.amp_relerr = ", demo_result.amp_relerr)
println("demo_result.pattern_relerr = ", demo_result.pattern_relerr)
