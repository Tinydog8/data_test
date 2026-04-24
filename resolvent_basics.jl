#!/usr/bin/env julia

#=
最小 Julia 学习脚本：按代码块逐步运行的 resolvent analysis 入门
===============================================================

这份脚本专门改成“不要封装成函数，而是分块顺序执行”的形式。
如果你使用 VS Code 的 Julia 插件，可以直接按 `# %%` 单元逐块运行。

核心目标只有一个：

    让你先完全看懂：
    谐波强迫 -> resolvent 算子 -> SVD -> 最优 forcing / response

依旧只保留最基础模型：

    ∂u/∂t = ν ∂²u/∂y² + d(y, t),    y ∈ (0, 1)
    u(0, t) = 0,  u(1, t) = 0

如果令

    d(y, t) = d̂(y; ω) exp(-iωt),
    u(y, t) = û(y; ω) exp(-iωt),

那么可得频域方程

    (iωI - ν∂yy) û = d̂

于是 resolvent operator 为

    H(ω) = (iωI - ν∂yy)^(-1)

在满足边界条件的标准化正弦基

    ϕ_n(y) = √2 sin(nπy)

下，有

    ∂yy ϕ_n = -(nπ)^2 ϕ_n

因此第 n 个模态对应的 resolvent 增益为

    H_n(ω) = 1 / (iω + ν(nπ)^2)
    σ_n(ω) = |H_n(ω)|
           = 1 / sqrt(ω^2 + [ν(nπ)^2]^2)

这正是最适合入门和验证的极简案例。
=#

# %%
# --------------------------------------------------------------------
# 第 0 块：加载库
# --------------------------------------------------------------------
# 这里只保留最必要的两个标准库：
# - LinearAlgebra: 用于 SVD 和范数
# - Printf: 用于更整齐地打印结果

using LinearAlgebra
using Printf


# %%
# --------------------------------------------------------------------
# 第 1 块：设置最少参数
# --------------------------------------------------------------------
# 这一块只做“实验设定”，建议先单独运行并观察变量。
#
# 你可以先改：
#   - nu    : 扩散/耗散强度
#   - omega : 谐波强迫角频率
#   - nmodes: 保留多少个正弦模态
#
# 这个模型里没有网格离散误差，因为我们直接在正弦本征基里工作。
# nmodes 只是“保留多少个模态做展示”，不是 PDE 离散精度的来源。

nu     = 0.02
omega  = 1.00
nmodes = 8

# 物理空间网格仅用于后面把模态重构成曲线形状，不参与求解。
ygrid = range(0.0, 1.0; length = 201)

println("参数设置完成：")
@printf("nu     = %.6f\n", nu)
@printf("omega  = %.6f\n", omega)
@printf("nmodes = %d\n", nmodes)


# %%
# --------------------------------------------------------------------
# 第 2 块：写出 Laplacian 在正弦基下的本征值
# --------------------------------------------------------------------
# 边界条件 u(0)=u(1)=0 对应最自然的基底：
#
#     ϕ_n(y) = √2 sin(nπy),   n = 1,2,3,...
#
# 且
#
#     ∂yy ϕ_n = -(nπ)^2 ϕ_n
#
# 因此 Laplacian 的“正特征值大小”记作
#
#     λ_n = (nπ)^2
#
# 在这个基底下，原方程完全解耦成若干个标量代数方程。

lambdas = [(n * π)^2 for n in 1:nmodes]

println("\n每个模态对应的 λ_n = (nπ)^2：")
for n in 1:nmodes
    @printf("mode %d : lambda = %.10f\n", n, lambdas[n])
end


# %%
# --------------------------------------------------------------------
# 第 3 块：构造 resolvent 矩阵 H(ω)
# --------------------------------------------------------------------
# 频域问题：
#
#     (iωI - ν∂yy) û = d̂
#
# 由于
#
#     ∂yy ϕ_n = -(nπ)^2 ϕ_n = -λ_n ϕ_n
#
# 所以第 n 个模态上有
#
#     (iω + νλ_n) û_n = d̂_n
#
# 因而
#
#     û_n = H_n d̂_n
#     H_n = 1 / (iω + νλ_n)
#
# 也就是说，resolvent H 在这组基底下是一个对角矩阵。

diagonal_entries = ComplexF64[]
for λ in lambdas
    push!(diagonal_entries, 1 / (im * omega + nu * λ))
end

H = Diagonal(diagonal_entries)

println("\nresolvent 对角元 H_n：")
for n in 1:nmodes
    @printf("mode %d : Re(H_n) = % .10e, Im(H_n) = % .10e\n",
            n, real(H[n, n]), imag(H[n, n]))
end


# %%
# --------------------------------------------------------------------
# 第 4 块：写出解析增益，并与 resolvent 的 SVD 对照
# --------------------------------------------------------------------
# 由于 H 是对角矩阵，且各模态独立，
# 第 n 个模态的 singular value 就是 |H_n|：
#
#     σ_n = |H_n|
#         = 1 / sqrt(ω^2 + (νλ_n)^2)
#
# 这是我们最重要的解析 benchmark。

exact_gains = Float64[]
for λ in lambdas
    push!(exact_gains, 1 / sqrt(omega^2 + (nu * λ)^2))
end

# 数值上直接对 H 做 SVD。
#
# H = U Σ V*
#
# 含义：
# - V 的列向量：forcing 空间中的最优方向
# - U 的列向量：response 空间中的最优方向
# - Σ 的对角元：各方向上的放大率

F = svd(Matrix(H))
numeric_gains = F.S

println("\n模态增益对照（数值 SVD vs 解析公式）")
println("--------------------------------------------------------------")
@printf("%6s  %18s  %18s  %14s\n", "mode", "numeric sigma", "exact sigma", "rel. error")
for n in 1:nmodes
    relerr = abs(numeric_gains[n] - exact_gains[n]) / exact_gains[n]
    @printf("%6d  %18.10e  %18.10e  %14.6e\n",
            n, numeric_gains[n], exact_gains[n], relerr)
end


# %%
# --------------------------------------------------------------------
# 第 5 块：找出“最优 forcing”与“最优 response”
# --------------------------------------------------------------------
# 对本问题而言，理论上最大增益必然出现在最低模态 n=1。
# 所以：
#
#     最优 forcing  = 第一正弦模态
#     最优 response = H 乘以该 forcing 后得到的响应
#
# 这里先在“模态系数空间”中表示：
#
#     forcing_coeffs  = [1, 0, 0, ..., 0]^T
#     response_coeffs = H * forcing_coeffs
#
# 这一步最适合帮助理解：
# “forcing 的方向”和“response 的方向”在基础问题里其实非常清楚。

optimal_forcing_coeffs = zeros(ComplexF64, nmodes)
optimal_forcing_coeffs[1] = 1.0 + 0.0im

optimal_response_coeffs = H * optimal_forcing_coeffs

gain_from_direct_solve = norm(optimal_response_coeffs) / norm(optimal_forcing_coeffs)

println("\n最重要的三个量：")
@printf("最大奇异值 sigma_1      = %.10e\n", numeric_gains[1])
@printf("直接计算得到的响应增益   = %.10e\n", gain_from_direct_solve)
@printf("二者之差                 = %.6e\n",
        abs(numeric_gains[1] - gain_from_direct_solve))


# %%
# --------------------------------------------------------------------
# 第 6 块：把模态系数重构到物理空间
# --------------------------------------------------------------------
# 前面我们在“模态空间”工作，现在把 forcing / response 还原回物理空间：
#
#     f(y) = Σ a_n ϕ_n(y)
#     u(y) = Σ b_n ϕ_n(y)
#
# 其中
#
#     ϕ_n(y) = √2 sin(nπy)
#
# 这一块故意不用函数封装，而是直接展开写循环，方便你逐行看懂。

forcing_profile  = zeros(ComplexF64, length(ygrid))
response_profile = zeros(ComplexF64, length(ygrid))

for j in eachindex(ygrid)
    y = ygrid[j]

    forcing_sum  = 0.0 + 0.0im
    response_sum = 0.0 + 0.0im

    for n in 1:nmodes
        phi_n = sqrt(2.0) * sin(n * π * y)
        forcing_sum  += optimal_forcing_coeffs[n]  * phi_n
        response_sum += optimal_response_coeffs[n] * phi_n
    end

    forcing_profile[j]  = forcing_sum
    response_profile[j] = response_sum
end

println("\n已完成 forcing / response 的物理空间重构。")


# %%
# --------------------------------------------------------------------
# 第 7 块：打印少量采样点，观察最优形状
# --------------------------------------------------------------------
# 理论上：
# - 最优 forcing 应该就是第一正弦模态 sin(πy)
# - response 与它形状相同，只是乘上一个复数增益 H_1
#
# 这里不做绘图，只打印若干点，保持脚本极简。

println("\n物理空间采样值（实部）")
println("--------------------------------------------------------------")
@printf("%8s  %18s  %18s\n", "y", "forcing(real)", "response(real)")

sample_ids = [1, 41, 81, 101, 121, 161, 201]
for idx in sample_ids
    @printf("%8.3f  %18.10e  %18.10e\n",
            ygrid[idx], real(forcing_profile[idx]), real(response_profile[idx]))
end


# %%
# --------------------------------------------------------------------
# 第 8 块：结果解释
# --------------------------------------------------------------------
# 这一块不做新计算，只把前面的结果整理成可以直接学习的结论。

println("\n结论说明")
println("--------------------------------------------------------------")
println("1. 这个最小模型中的 resolvent 算子在正弦基下完全对角化。")
println("2. 因此每个模态彼此独立，singular value 可直接由解析公式得到。")
println("3. 最优 forcing 是第一模态，因为它的耗散最弱。")
println("4. 当 omega 增大时，所有模态响应都会减弱。")
println("5. 当 nu 增大时，高阶模态会更快被抑制。")
println()
println("把它升级到论文问题时，需要逐步加入：")
println("  - 状态向量而不是标量变量")
println("  - 平均流、波数和耦合项")
println("  - 输入矩阵 B、输出矩阵 C")
println("  - 对每个 (kx, kz, omega) 的 resolvent 做 SVD")
println()
println("这份脚本的目的不是复现论文，而是让你逐块跑通最基础的 resolvent 思想。")
