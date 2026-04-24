#!/usr/bin/env julia

#=
最小 Julia 学习脚本：1D 受迫扩散方程的 resolvent analysis
===========================================================

这份脚本故意只保留论文方法里“最基础、最核心”的部分：

    线性系统 + 持续谐波强迫 + 频率响应 + 最大放大

而把下列复杂内容全部拿掉：

    - 三维流动
    - Orr-Sommerfeld / Squire 耦合
    - 波流相互作用中的 Stokes drift
    - 涡黏性剖面
    - LES 数据拟合
    - 绘图、文件输出、命令行解析等无关功能

----------------------------------------------------------------------
一、为什么用这个模型？
----------------------------------------------------------------------

论文研究的是更复杂的 resolvent 问题：

    û = T(kx, kz, ω) d̂

其中 T 是一个由线性化 Navier-Stokes / CL 方程构成的传递算子。

要学会它，最关键的不是一上来就处理复杂流动，而是先真正理解：

    1. 什么是“谐波强迫”
    2. 什么是“频率响应 / resolvent operator”
    3. 为什么要对这个算子做 SVD
    4. 什么叫“最容易被放大的 forcing / response”

下面的极简模型已经完整包含了这四点。

----------------------------------------------------------------------
二、所求解的基础方程
----------------------------------------------------------------------

考虑区间 y ∈ (0, 1) 上的一维受迫扩散方程：

    ∂u/∂t = ν ∂²u/∂y² + d(y, t),

配 Dirichlet 边界条件：

    u(0, t) = 0,   u(1, t) = 0.

这里：

    u(y, t) : 系统响应
    d(y, t) : 外部强迫
    ν       : 扩散系数（对应流体问题里的黏性/耗散）

如果施加单频谐波强迫：

    d(y, t) = d̂(y; ω) exp(-i ω t),

并设系统响应也具有相同频率：

    u(y, t) = û(y; ω) exp(-i ω t),

代回原方程后得到频域问题：

    (i ω I - ν ∂yy) û = d̂.

于是 resolvent operator 定义为：

    H(ω) = (i ω I - ν ∂yy)^(-1),

所以有：

    û = H(ω) d̂.

这就是论文中 resolvent 公式最本质的原型。

----------------------------------------------------------------------
三、为什么这个模型还能“精确验证”？
----------------------------------------------------------------------

由于边界条件是 u(0)=u(1)=0，最自然的基底是正弦基：

    ϕ_n(y) = √2 sin(n π y),   n = 1, 2, 3, ...

它满足：

    ∂yy ϕ_n = -(n π)^2 ϕ_n.

因此在这个基底下，算子完全对角化：

    H_n(ω) = 1 / (i ω + ν (n π)^2).

于是第 n 个模态的放大量（也就是 singular value）有解析表达式：

    σ_n(ω) = |H_n(ω)|
           = 1 / sqrt(ω^2 + [ν (n π)^2]^2).

这意味着：

    - 最大放大一定来自 n = 1 的最低模态
    - 数值结果可以直接对照解析结果
    - 这是学习 resolvent 最理想的基础例子

----------------------------------------------------------------------
四、这和论文的关系是什么？
----------------------------------------------------------------------

这份代码是“论文方法的最小骨架”：

    当前脚本:
        标量 PDE  ->  一个 resolvent 算子 H(ω)

    论文中的系统:
        向量 PDE  ->  一个矩阵微分算子 T(kx, kz, ω)

真正升级到论文时，你只需要逐步增加三类复杂度：

    (1) 把标量 u 改成状态向量 q = [v, η] 或速度/涡量组合
    (2) 把 ν ∂yy 改成含平均流、波数、Stokes drift 的线性算子
    (3) 把单纯 H(ω) 改成输入-输出形式 T = C(iωE - F)^(-1)B

但“谐波强迫 -> resolvent -> SVD -> 最优响应”的核心逻辑完全不变。

=#

using LinearAlgebra
using Printf

# --------------------------------------------------------------------
# 1. 构造正弦基下的最小 resolvent 模型
# --------------------------------------------------------------------

"""
    build_resolvent(nu, omega, nmodes)

在正弦基 ϕ_n(y)=√2 sin(nπy) 下构造 resolvent 矩阵。

数学上：

    H_n(ω) = 1 / (iω + ν(nπ)^2)

因为不同正弦模态彼此独立，所以 H 在该基底下是对角矩阵。

返回值：

    H        : resolvent 矩阵（复数对角矩阵）
    lambdas  : 每个模态对应的 Laplacian 本征值 (nπ)^2
"""
function build_resolvent(nu::Float64, omega::Float64, nmodes::Int)
    lambdas = [(n * π)^2 for n in 1:nmodes]

    # 对角元素对应:
    #     H_n = 1 / (iω + ν λ_n)
    # 其中 λ_n = (nπ)^2
    diagonal_entries = ComplexF64[
        1 / (im * omega + nu * λ) for λ in lambdas
    ]

    H = Diagonal(diagonal_entries)
    return H, lambdas
end


"""
    analytic_gains(nu, omega, lambdas)

给出每个模态增益的解析表达式：

    σ_n(ω) = 1 / sqrt(ω^2 + [ν λ_n]^2)
"""
function analytic_gains(nu::Float64, omega::Float64, lambdas::Vector{Float64})
    return [1 / sqrt(omega^2 + (nu * λ)^2) for λ in lambdas]
end


"""
    sine_mode(n, y)

标准化正弦基函数：

    ϕ_n(y) = √2 sin(nπy)

选这个标准化是为了使 L2 内积下各模态正交归一，便于解释
“forcing 的单位能量”和“response 的放大量”。
"""
@inline function sine_mode(n::Int, y::Float64)
    return sqrt(2.0) * sin(n * π * y)
end


"""
    reconstruct_field(coeffs, ygrid)

将模态系数重构为物理空间中的函数值：

    u(y) = Σ coeffs[n] ϕ_n(y)

这里 coeffs 可以是 forcing 的模态系数，也可以是 response 的模态系数。
"""
function reconstruct_field(coeffs::AbstractVector{ComplexF64}, ygrid)
    values = zeros(ComplexF64, length(ygrid))
    for (j, y) in enumerate(ygrid)
        acc = 0.0 + 0.0im
        for n in eachindex(coeffs)
            acc += coeffs[n] * sine_mode(n, y)
        end
        values[j] = acc
    end
    return values
end


# --------------------------------------------------------------------
# 2. 主程序：构造模型、做 SVD、对照解析结果
# --------------------------------------------------------------------

function main()
    # -----------------------------
    # 可直接修改的最少参数
    # -----------------------------
    nu     = 0.02     # 扩散系数/耗散强度
    omega  = 1.00     # 谐波强迫角频率
    nmodes = 8        # 截断到前 nmodes 个正弦模态

    # 物理空间网格只用于把模态系数重构成函数值，便于理解结果。
    # 它不参与求解，因此不会影响数值精度。
    ygrid = range(0.0, 1.0; length = 201)

    # 构造 resolvent
    H, lambdas = build_resolvent(nu, omega, nmodes)

    # 对 resolvent 做 SVD：
    #     H = U Σ V*
    #
    # 其中：
    #     V 的第一列 = 最优 forcing 方向
    #     U 的第一列 = 最优 response 方向
    #     Σ[1]       = 最大放大量
    #
    # 对本问题而言，H 是对角矩阵，因此最优 forcing 就是第一模态。
    F = svd(Matrix(H))
    numeric_gains = F.S
    exact_gains = analytic_gains(nu, omega, lambdas)

    # 由于本问题的最优 forcing 就是第一正弦模态，下面显式构造它，
    # 以便把“模态空间中的向量”重构到物理空间中。
    optimal_forcing_coeffs = zeros(ComplexF64, nmodes)
    optimal_forcing_coeffs[1] = 1.0 + 0.0im

    # 对应的最优响应为:
    #     û = H f̂
    optimal_response_coeffs = H * optimal_forcing_coeffs

    # 数值验证：
    # 对单位范数 forcing，响应范数应等于最大 singular value。
    gain_from_direct_solve = norm(optimal_response_coeffs) / norm(optimal_forcing_coeffs)

    # 重构物理空间中的 forcing / response 形状。
    forcing_profile  = reconstruct_field(optimal_forcing_coeffs, ygrid)
    response_profile = reconstruct_field(optimal_response_coeffs, ygrid)

    # -----------------------------
    # 输出结果
    # -----------------------------
    println("==============================================================")
    println("最小 resolvent 学习示例：1D 受迫扩散方程")
    println("==============================================================")
    @printf("nu     = %.6f\n", nu)
    @printf("omega  = %.6f\n", omega)
    @printf("nmodes = %d\n", nmodes)
    println()

    println("模态增益对照（数值 SVD vs 解析公式）")
    println("--------------------------------------------------------------")
    @printf("%6s  %18s  %18s  %14s\n", "mode", "numeric sigma", "exact sigma", "rel. error")
    for n in 1:nmodes
        relerr = abs(numeric_gains[n] - exact_gains[n]) / exact_gains[n]
        @printf("%6d  %18.10e  %18.10e  %14.6e\n",
                n, numeric_gains[n], exact_gains[n], relerr)
    end
    println()

    println("最重要的三个结论")
    println("--------------------------------------------------------------")
    @printf("1) 最大放大量 sigma_1 = %.10e\n", numeric_gains[1])
    @printf("2) 直接求解得到的增益    = %.10e\n", gain_from_direct_solve)
    @printf("3) 二者之差              = %.6e\n",
            abs(numeric_gains[1] - gain_from_direct_solve))
    println()

    println("物理解释")
    println("--------------------------------------------------------------")
    println("1. 最优 forcing 是第一正弦模态 sin(pi y)。")
    println("2. 原因是低阶模态的扩散衰减最弱，因此最容易被持续谐波强迫放大。")
    println("3. 当 omega 变大时，所有模态的增益都会下降。")
    println("4. 当 nu 变大时，高阶模态会被更强烈地抑制。")
    println()

    println("给出少量物理空间采样值，帮助理解“最优 forcing / response 形状”")
    println("--------------------------------------------------------------")
    @printf("%8s  %18s  %18s\n", "y", "forcing(real)", "response(real)")
    sample_ids = [1, 41, 81, 101, 121, 161, 201]
    for idx in sample_ids
        @printf("%8.3f  %18.10e  %18.10e\n",
                ygrid[idx], real(forcing_profile[idx]), real(response_profile[idx]))
    end
    println()

    println("如何把这份基础代码升级到论文问题？")
    println("--------------------------------------------------------------")
    println("步骤 1: 把标量 u 改成状态向量 q（例如速度/涡量变量）。")
    println("步骤 2: 把标量算子 nu*d²/dy² 改成含平均流和波数的线性算子。")
    println("步骤 3: 把 H(omega) 改成 T(kx, kz, omega) 的输入-输出形式。")
    println("步骤 4: 再对每个 (kx, kz, omega) 做 SVD，寻找主导放大结构。")
    println()

    println("这份脚本的用途不是复现论文，而是让你先把 resolvent 的最小骨架学扎实。")
end


main()
