## 公式—代码对应说明（闭环传递函数与 H∞ 范数）

本文档把你推导的频域方程，与我在 `closed_loop_transfer.jl` 中给出的“块矩阵一次求解”实现逐一对应起来，重点解释传递函数 **G** 在代码中到底是什么。

---

### 1. 连续方程（频域）与符号

我们考虑展向傅里叶模态（流向常数，$k_x=0$），对时间做拉普拉斯/傅里叶变换，令

- $s=\sigma+i\omega$（做稳态频响时取 $s=i\omega$）
- $D = \frac{\partial}{\partial z}$
- $\Delta = D^2 - k^2$

在你推导的 OS–Squire 形式（包含 vortex force 的耦合）下，频域线性系统（保留强迫）可写为：

\[
\begin{aligned}
&(sI-\Delta/Re)\,\hat u + U'(z)\,\hat w = \hat f_x,\\
&-\frac{1}{La}k^2\,u_s'(z)\,\hat u + (s\Delta-\Delta^2/Re)\,\hat w
= (-ikD)\hat f_y + (-k^2 I)\hat f_z.
\end{aligned}
\]

其中 $u_s'(z)=Du_s(z)$。

---

### 2. 块矩阵形式：把“闭环耦合”写成一次线性求解

将未知量与输入堆叠为

\[
\psi = \begin{bmatrix}\hat u\\ \hat w\end{bmatrix},\qquad
\hat{\mathbf f} = \begin{bmatrix}\hat f_x\\ \hat f_y\\ \hat f_z\end{bmatrix}.
\]

则上面的两条方程可写成

\[
\underbrace{
\begin{bmatrix}
 sI-\Delta/Re & U'\\
 -(1/La)k^2\,\mathrm{diag}(u_s') & s\Delta-\Delta^2/Re
\end{bmatrix}}_{A(s)}
\underbrace{\begin{bmatrix}\hat u\\ \hat w\end{bmatrix}}_{\psi}
=
\underbrace{
\begin{bmatrix}
 I & 0 & 0\\
 0 & -ikD & -k^2 I
\end{bmatrix}}_{B}
\underbrace{\begin{bmatrix}\hat f_x\\ \hat f_y\\ \hat f_z\end{bmatrix}}_{\hat{\mathbf f}}.
\]

**关键点**：你之前用 $G_1,G_2,G_3,\Gamma$ 表达“闭环”，在离散与边界条件（tau 行替换）存在时，很容易出现算子次序/BC 不一致导致的不等价。块矩阵一次求解等价于“先消元再闭环”，但数值上更稳、更不易写错。

---

### 3. 传递函数 \(G(s)\) 的严格定义（本文重点）

我们关心输出是壁法向速度 $\hat w$，输入是 $\hat{\mathbf f}$。由块矩阵方程

\[
\psi = A(s)^{-1}B\,\hat{\mathbf f}.
\]

取 $\psi$ 的下半块就是 $\hat w$。令输出选择矩阵

\[
C = \begin{bmatrix}0 & I\end{bmatrix},
\]

则

\[
\hat w = C\,A(s)^{-1}B\,\hat{\mathbf f} \equiv G(s)\,\hat{\mathbf f}.
\]

因此

\[
\boxed{\;G(s) = C\,A(s)^{-1}B\;}
\]

在代码里，**返回值就是离散后的矩阵版 $G(i\omega)$**，维度为 $n\times 3n$：

- 行：$w$ 在 $n$ 个 Chebyshev 网格点上的值
- 列：堆叠输入 $[f_x;f_y;f_z]$ 的 $3n$ 个自由度

---

### 4. 代码对应关系（逐块对照）

下面对照 `closed_loop_transfer.jl` 的实现。

#### 4.1 构造 \(\Delta\) 与 \(\Delta^2\)

公式：
\[
\Delta = D^2 - k^2I,\qquad \Delta^2 = \Delta\,\Delta.
\]

代码（`build_matrices`）：

```julia
D2 = D1 * D1
Δ  = D2 .- k^2 .* I
Δ2 = Δ * Δ
```

> 重要修正：`z` 反转时必须用同一个置换同步重排 `D1`，否则导数方向会错。

#### 4.2 离散 \(u_s'(z)\)（必须先求导再对角乘）

公式：
\[
\text{耦合项} \propto \mathrm{diag}(u_s'(z))\,\hat u.
\]

代码（`G_closed_frequency`）：

```julia
u_s′ = D1 * u_s_vec
Us′diag = Diagonal(ComplexF64.(u_s′))
```

这对应的是 $u_s'(z)=Du_s(z)$，再形成对角矩阵乘以 $u$。

#### 4.3 构造块算子 \(A(s)\)

公式：
\[
A(s)=\begin{bmatrix}
 sI-\Delta/Re & U'\\
 -(1/La)k^2\,\mathrm{diag}(u_s') & s\Delta-\Delta^2/Re
\end{bmatrix}.
\]

代码：

```julia
Auu = s .* I .- Δc ./ Re
Auw = U′diag
Awu = -(1.0 / La) * k^2 .* Us′diag
Aww = s .* Δc .- Δ2c ./ Re
A = [Auu  Auw;
     Awu  Aww]
```

#### 4.4 构造输入算子 \(B\)

公式：
\[
B = \begin{bmatrix} I & 0 & 0\\ 0 & -ikD & -k^2I\end{bmatrix}.
\]

代码：

```julia
Fx = I
Fy = (-im * k) .* D1c
Fz = (-k^2) .* I
RHS = [Fx  Z   Z;
       Z   Fy  Fz]
```

其中 `RHS` 就是离散后的 $B$。

#### 4.5 边界条件（tau 行替换）

你使用的边界条件为

- $u(\pm 1)=0$
- $w(\pm 1)=0$
- $Dw(\pm 1)=0$

代码通过**替换相应方程行**为边界约束，并把对应 RHS 置零实现。例如：

- `u(-1)=0`：把 $u$ 方程的第 1 行替换为 $u_1=0$。
- `w(-1)=0`：把 $w$ 方程的第一行（在块矩阵里是第 `n+1` 行）替换为 $w_1=0$。
- `Dw(-1)=0`：把第 `n+2` 行替换为 $D1[1,:] * w = 0$。

（完整实现见 `G_closed_frequency` 的 “Enforce BC” 段。）

#### 4.6 求解并提取 \(G(i\omega)\)

由定义：
\[
G(s)=C\,A(s)^{-1}B.
\]

代码：

```julia
F = lu(A) \ RHS          # F ≈ A(s)^{-1} B
Wmap = F[n+1:2n, :]      # 取下半块：w
return Wmap              # 这就是 G(iω)
```

这里的 `Wmap` 就是本文定义的 **传递函数矩阵 $G(i\omega)$**。

---

### 5. 无穷范数 \(\|G\|_\infty\) 在代码中对应什么

公式（离散频率扫描近似）：
\[
\|G\|_\infty \approx \max_{\omega\in\Omega}\ \sigma_{\max}(G(i\omega)).
\]

代码（`hinf_norm`）：

```julia
σ = maximum(svdvals(Matrix(G)))
maxσ = max(maxσ, σ)
```

- `svdvals(G)` 产生所有奇异值
- `maximum(...)` 取 $\sigma_{\max}$
- 扫描 `ω_range` 取最大值即近似 $\|G\|_\infty$

---

### 6. 和你原来的 \(G_1,G_2,G_3,\Gamma\) 写法的关系（简述）

你推导的“闭环形式”本质是对块系统做 Schur 消元得到：

\[
\hat w = (I+\Gamma)^{-1}\big(G_3\hat f_x + G_1\hat f_y + G_2\hat f_z\big),
\]

其中 $G_1,G_2,G_3,\Gamma$ 都是算子（矩阵）且乘法次序不可随意交换。

本文代码等价地直接计算 $C A^{-1} B$，避免了显式构造 $(I+\Gamma)^{-1}$ 时的次序与边界条件细节问题。

---

### 7. 最小使用示例

```julia
include("closed_loop_transfer.jl")

N = 200
Re = 1.0
La = 0.1
flow_type = "couette"

k_values = range(0.1, 9.0, length=90)
ω_range = range(-5.0, 5.0, length=201)

k_vals, norms, bestωs = compute_hinf_vs_k(
    k_values; flow_type=flow_type, N=N, Re=Re, La=La,
    u_s_params=(1.0, 2pi/2.4), ω_range=ω_range
)

println(maximum(norms), " at k=", k_vals[argmax(norms)], " ω=", bestωs[argmax(norms)])
```
