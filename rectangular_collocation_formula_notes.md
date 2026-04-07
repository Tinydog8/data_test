# `rectangular_collocation` 代码与论文公式对应说明

本文档说明你给出的 `rectangular_collocation` 实现与论文  

> Xuan & Shen, *Resolvent model-based analyses of coherent structures in Langmuir turbulence*, JFM (2025)

中的公式、附录 B 离散方法之间的对应关系。

文中最核心的连续方程和算子定义位于：

- 正文 §2.1，公式 (2.5)–(2.19)
- 附录 B，公式 (B1)–(B6)

---

## 1. 连续问题：代码想离散的到底是什么

论文将状态变量定义为

\[
\xi = \begin{bmatrix} v \\ \omega_y \end{bmatrix},
\]

其中：

- \(v\)：垂向速度扰动
- \(\omega_y\)：垂向涡量扰动

在对 \((x,z,t)\) 做 Fourier 变换后，单个 \((k_x,k_z,\omega)\) 模式满足论文公式 (2.9)：

\[
-(i\omega E + F)\hat{\xi} = B \hat d .
\]

并通过论文公式 (2.17)

\[
\hat u = C \hat \xi
\]

得到三分量速度响应 \(\hat u = [\hat u,\hat v,\hat w]^T\)。

文中也写成传递算子形式（论文公式 (2.19)）：

\[
\hat u = T \hat d, \qquad
T = C(i\omega E - F)^{-1}B .
\]

你的代码目标就是离散这个 \(T\)。

---

## 2. 状态空间与论文中的 \(v,\omega_y\) 变量

### 代码中的相关对象

- `v` 对应论文中的 \(v\)
- `ω`（代码里实际用 `w` 命名的那一组未知量）对应论文中的 \(\omega_y\)

注意：

- 代码里的 `y_w`, `Dw`, `D2w` **不是 spanwise velocity \(w\) 的网格**
- 它们实际上对应的是 **\(\omega_y\) 方程的网格**

因此代码中：

- `Uv`, `Usv`, `nuTv`：表示定义在 \(v\)-equation 网格上的剖面
- `Uw`, `Usw`, `nuTw`：表示定义在 \(\omega_y\)-equation 网格上的同一个物理剖面

它们不是不同物理分量，而是同一函数在两套离散网格上的采样。

---

## 3. `RectGrid`：附录 B 的矩形谱配置离散

### 3.1 论文对应

附录 B 说明：

- 对一个 \(m\) 阶方程，未知函数离散在 \(N+m\) 个 second-kind Chebyshev 点上
- PDE 本体在 \(N\) 个 first-kind Chebyshev 点上 enforce
- 再附加 \(m\) 个边界条件

对应附录 B 的公式 (B3) 与 (B4)：

\[
x_j = -\cos\left(\frac{j\pi}{N+m-1}\right), \qquad j=0,\dots,N+m-1
\]

\[
\check x_j = -\cos\left(\frac{(j+1/2)\pi}{N}\right), \qquad j=0,\dots,N-1
\]

---

### 3.2 `RectGrid` 中的字段对应

```julia
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
```

对应关系如下。

#### `N`
正文和附录 B 中 PDE 被 enforce 的第一类 Chebyshev 内点个数。

#### `Nv = N + 4`
因为 \(v\)-equation 是四阶。

#### `Nw = N + 2`
因为 \(\omega_y\)-equation 是二阶。

#### `xi_v`, `xi_w`
second-kind Chebyshev 点，对应附录 B 的 (B3)。

#### `xi_int`
first-kind Chebyshev 点，对应附录 B 的 (B4)。

#### `y_v`, `y_w`, `y_int`
把参考区间 \([-1,1]\) 映射到物理区间 \([-H,0]\) 的垂向坐标。

这里代码使用

\[
y = -\frac{H}{2}(x+1)
\]

所以：

- \(x=1 \mapsto y=-H\)
- \(x=-1 \mapsto y=0\)

**注意：这意味着当前实现的数组顺序是“海底 \(\to\) 海表”。**

这点与你有些脚本中的中文注释“表面在首索引”不一致，使用时一定要统一。

#### `Dv`, `D2v`, `Dw`, `D2w`
在 second-kind Chebyshev 点上的一、二阶微分矩阵。

#### `Pv`, `Pw`
从 second-kind 点重采样到 first-kind 内点的插值矩阵，对应附录 B 中的 \(P\)。

#### `I_vw`, `I_wv`
在两套 second-kind 点之间做插值：

- `I_vw`: \(\omega_y\)-grid \(\to\) \(v\)-grid
- `I_wv`: \(v\)-grid \(\to\) \(\omega_y\)-grid

这是为了在两个状态变量位于不同网格时构造耦合项。

#### `I_Nv_N`, `I_Nw_N`
将 first-kind interior forcing 网格插值到 second-kind 网格，用于组装输入算子 \(B\)。

#### `w_y_int`
离散的能量权重，用于近似论文公式 (2.21) 中的能量内积

\[
\langle f,g \rangle_E = \int_{-H}^0 g^* f\,dy .
\]

你当前实现使用 \(\theta\)-空间中点法则：

\[
x = -\cos\theta,\qquad
\int_{-1}^{1} f(x)\,dx
=
\int_0^\pi f(-\cos\theta)\sin\theta\,d\theta
\]

在

\[
\theta_j = \frac{(j-1/2)\pi}{N}
\]

上得到

\[
w_j^{(\xi)} \approx \frac{\pi}{N}\sin\theta_j,
\qquad
w_j^{(y)} = \frac{H}{2} w_j^{(\xi)} .
\]

---

## 4. `trefethen_cheb_dmat`：Chebyshev-Lobatto 微分矩阵

```julia
function trefethen_cheb_dmat(N::Int)
```

这个函数返回：

- Lobatto 节点 \(x_j = \cos(j\pi/(N-1))\)
- 一阶导数矩阵 \(D\)
- 二阶导数矩阵 \(D^2\)

这是经典 Trefethen Chebyshev collocation 公式。

在当前代码中，这些参考区间上的微分矩阵后续还会乘以坐标缩放因子

\[
\frac{d}{dy} = -\frac{2}{H}\frac{d}{dx}
\]

因此：

```julia
sc = -2 / H
Dv .*= sc
D2v .*= sc^2
```

是把参考导数矩阵变成物理坐标 \(y\) 下的导数矩阵。

---

## 5. `cheb1_points`：第一类 Chebyshev 内点

```julia
cheb1_points(N::Int) = [-cos(pi * (j + 0.5) / N) for j in 0:N-1]
```

这正是附录 B 的 (B4)：

\[
\check x_j = -\cos\left(\frac{(j+1/2)\pi}{N}\right).
\]

这些点是 PDE 本体被 enforce 的位置。

---

## 6. `barycentric_interp_matrix`：附录 B 中的重采样矩阵 \(P\)

附录 B 说明矩形离散中的核心对象是：

\[
P D^k
\]

其中：

- \(D^k\) 是 second-kind 点上的导数矩阵
- \(P\) 把 second-kind 上的多项式值重采样到 first-kind 内点

代码中的：

```julia
Pv = barycentric_interp_matrix(collect(xi_v), collect(xi_int))
Pw = barycentric_interp_matrix(collect(xi_w), collect(xi_int))
```

正是在构造这个 \(P\)。

此外：

```julia
I_vw = barycentric_interp_matrix(collect(xi_w), collect(xi_v))
I_wv = barycentric_interp_matrix(collect(xi_v), collect(xi_w))
```

是在两套 second-kind 网格之间重采样。

---

## 7. `build_L_blocks`：论文中的 \(L_{OS},L_{Sq}\) 与耦合项

论文给出的频域算子块为：

论文公式 (2.11) 为：
\[
F=
\begin{bmatrix}
L_{OS} & -ik_z U^{s\prime}\\
-ik_z U' & L_{Sq}
\end{bmatrix}.
\]

其中

论文公式 (2.13) 为：
\[
L_{OS}
=
-ik_x U^L \hat\Delta
+ik_x U''
+\nu_T \hat\Delta^2
+2\nu_T' D\hat\Delta
+\nu_T''(D^2+k^2I),
\]

论文公式 (2.14) 为：
\[
L_{Sq}
=
-ik_x U^L
+\nu_T \hat\Delta
+\nu_T' D.
\]

---

### 7.1 `Delta_hat`

```julia
Delta_hat(D2::AbstractMatrix, k::Real) = complex.(D2) - k^2 * I
```

对应论文里的

\[
\hat\Delta = D^2 - k^2 I.
\]

---

### 7.2 `build_L_blocks`

```julia
function build_L_blocks(...)
```

这一函数的作用是构造：

- `Δv`：\(v\)-grid 上的 \(\hat\Delta\)
- `L_OS`：Orr-Sommerfeld 型块
- `L_Sq`：Squire 型块
- `F12`：上右耦合块 \(-ik_z U^{s\prime}\)
- `F21`：下左耦合块 \(-ik_z U'\)

代码中：

```julia
ULv = Uv .+ Usv
ULw = Uw .+ Usw
```

对应论文中的

\[
U^L = U + U^s.
\]

---

### 7.3 `L_OS`

代码：

```julia
L_OS = (-im * kx) .* (D_ULv * Δv) .+ (im * kx) .* D_Uppv .+ (D_nuTv * Δ2v) .+
       (2 .* D_dnupv * Dv_c * Δv) .+ (D_d2nupv * (D2v_c .+ k^2 .* I_Nv))
```

对应论文公式 (2.13)。

说明：

- `D_ULv * Δv` 对应 \(U^L \hat\Delta\)
- `D_Uppv` 对应 \(U''\)
- `D_nuTv * Δ2v` 对应 \(\nu_T \hat\Delta^2\)
- `D_dnupv * Dv_c * Δv` 对应 \(2\nu_T' D\hat\Delta\)
- `D_d2nupv * (D2v_c + k^2 I)` 对应代码作者采用的
  \(\nu_T''(D^2+k^2I)\) 实现形式

注意：论文排版/OCR 在这一项上最容易混乱，代码这里采用的是其频域算子写法。

---

### 7.4 `L_Sq`

代码：

```julia
L_Sq = (-im * kx) .* D_ULw .+ (D_nuTw * Δw) .+ (D_dnupw * Dw_c)
```

对应论文公式 (2.14)：

\[
L_{Sq} = -ik_x U^L + \nu_T \hat\Delta + \nu_T' D.
\]

---

### 7.5 `F12`, `F21`

代码：

```julia
F12 = (-im * kz) .* (diagc(dUsv) * g.I_vw)
F21 = (-im * kz) .* (diagc(dUw) * g.I_wv)
```

对应论文公式 (2.11) 的两个耦合块：

\[
F_{12} = -ik_z U^{s\prime},
\qquad
F_{21} = -ik_z U'.
\]

这里额外乘以 `I_vw` / `I_wv`，是因为两个变量在不同的垂向网格上，需要先做网格间插值。

---

## 8. `build_M_B_C_rect`：离散版的 \((i\omega E - F)\)、\(B\)、\(C\)

这个函数是整个代码的核心。

它构造：

- `M`：离散版的 \((i\omega E - F)\) 并附加边界条件
- `Bmat`：离散输入算子
- `Cm`：离散输出算子

---

### 8.1 `M`：离散系统矩阵

连续方程写成

\[
(i\omega E - F)\hat\xi = -B \hat d
\]

或等价地把符号吸收进右端。代码采用的是“把 interior PDE rows 和 BC rows 拼成一个方阵”的做法：

```julia
M[1:N, 1:Nv] = Pv * (im * omega .* Δv .- L_OS)
M[1:N, Nv+1:end] = -Pv * F12

M[Nv+1:Nv+N, 1:Nv] = -Pw * F21
M[Nv+1:Nv+N, Nv+1:end] = Pw * (im * omega .* I_Nw .- L_Sq)
```

其中：

- 前 \(N\) 行：\(v\)-equation 在 first-kind 内点上 enforce
- 再后 \(N\) 行：\(\omega_y\)-equation 在 first-kind 内点上 enforce
- 剩下 \(4+2\) 行：边界条件

这正是附录 B 的矩形系统思想：

\[
\begin{bmatrix}
P(\cdots) \\
L
\end{bmatrix}
u
=
\begin{bmatrix}
Pf\\
g
\end{bmatrix}.
\]

---

## 9. 边界条件在代码中的对应

论文原始 stress-free 边界条件在 \((v,\omega_y)\) 变量里是：

\[
\hat v = D^2\hat v = D\hat\omega_y = 0
\quad \text{at } y=0,-H.
\]

而你当前这版实现使用的是**混合边界条件**：

- 底部：no-slip
- 顶部：stress-free

代码中：

```julia
M[N+1, :] .= vbc([one(ComplexF64); zeros(ComplexF64, Nv - 1)])
M[N+2, :] .= vbc(D2v[1, :])
M[N+3, :] .= vbc(D2v[Nv, :])
M[N+4, :] .= vbc([zeros(ComplexF64, Nv - 1); one(ComplexF64)])
...
M[Nv+N+1, :] .= wbc(Dw[1, :])
M[Nv+Nw, :] .= wbc(Dw[Nw, :])
```

这版实现所对应的边界约束，需要根据数组端点与物理边界的对应仔细判读。

由于当前网格映射采用

\[
y = -\frac{H}{2}(x+1),
\]

且 `trefethen_cheb_dmat` 返回的 `x` 顺序是从 \(1\) 到 \(-1\)，因此：

- 数组首点对应海底 \(y=-H\)
- 数组末点对应海表 \(y=0\)

所以在解释这一段时，务必检查你最终使用的 `rectangular_collocation_noslip.jl` 与脚本注释是否一致。

---

## 10. `Bmat`：论文中的输入算子 \(B\)

论文公式 (2.12)：

\[
B=
\begin{bmatrix}
-ik_x D & -k^2 & -ik_z D\\
ik_z & 0 & -ik_x
\end{bmatrix}.
\]

你的代码中，先在 second-kind 网格上构造，再插值到 interior 点：

```julia
Binner = hcat((-im * kx) .* Dv, (-k^2) .* I_Nv, (-im * kz) .* Dv)
Bv = Pv * Binner * kron(Matrix{ComplexF64}(I, 3, 3), g.I_Nv_N)
Bw = Pw * hcat(im * kz .* I_Nw, Znw, (-im * kx) .* I_Nw) * kron(Matrix{ComplexF64}(I, 3, 3), g.I_Nw_N)
```

对应理解：

- `Bv`：进入 \(v\)-equation 的 forcing 块
- `Bw`：进入 \(\omega_y\)-equation 的 forcing 块

最后：

```julia
Bmat[1:N, :] .= Bv
Bmat[Nv+1:Nv+N, :] .= Bw
```

把两部分放到总系统右端。

---

## 11. `Cm`：论文中的输出算子 \(C\)

论文公式 (2.18)：

\[
C=\frac{1}{k^2}
\begin{bmatrix}
ik_x D & -ik_z \\
k^2 & 0 \\
ik_z D & ik_x
\end{bmatrix}.
\]

你的代码中：

```julia
row_u_v = (im * kx / k^2) .* (PvC * Dv)
row_u_w = (-im * kz / k^2) .* (PvC * Ivw)
row_v_v = PvC
row_w_v = (im * kz / k^2) .* (PvC * Dv)
row_w_w = (im * kx / k^2) .* (PvC * Ivw)
Cm = vcat(hcat(row_u_v, row_u_w), hcat(row_v_v, Znvw), hcat(row_w_v, row_w_w))
```

与论文公式完全对应，只是考虑了：

- \(v\) 与 \(\omega_y\) 位于不同网格
- 因此在 \(\omega_y\) 那一列上需要加 `Ivw` 做网格插值

---

## 12. `transfer_gain_rect`：离散传递算子与增益

```julia
function transfer_gain_rect(...)
    M, Bcopy, Cm = build_M_B_C_rect(...)
    X = M \ Bcopy
    T = Cm * X
    s1 = svdvals(T)[1]
    return (; G=abs2(s1), s1, T, M, Cm, g)
end
```

这正是在离散层面构造：

\[
T = C(i\omega E - F)^{-1} B
\]

并取最大奇异值：

\[
G = \sigma_1^2.
\]

对应论文：

\[
G(k_x,k_z,\omega)
=
\max_{\hat d \neq 0}
\frac{\|\hat u\|_E^2}{\|\hat d\|_E^2}
= \sigma_1^2.
\]

注意：`transfer_gain_rect` 返回的 `s1` 是 **Euclidean SVD** 下的最大奇异值。  
如果要与论文的能量内积严格一致，需要对 `T` 进行加权 SVD，即使用 `w_y_int`。

---

## 13. `weights_output_3N`：输出能量权重

```julia
function weights_output_3N(g::RectGrid)
    w = g.w_y_int
    vcat(w, w, w)
end
```

因为输出向量是

\[
\hat u = [\hat u,\hat v,\hat w]^T
\]

在离散后是长度 \(3N\) 的向量，所以离散能量内积的权重就是把同一套 \(y\)-积分权重复用三次。

这对应论文中的：

\[
\|f\|_E^2 = \int_{-H}^0 f^* f \,dy.
\]

---

## 14. 使用时最容易出错的几点

### 14.1 `Uv, Uw` 传入的必须是 Eulerian mean velocity \(U\)

代码注释已经强调：

```julia
# Uv,Uw = Eulerian mean U(y); Usv,Usw = Stokes drift. UL = U+Us (paper §2.2).
# Do not pass Lagrangian U^L as Uv or Stokes is double-counted.
```

也就是说：

- 如果你的原始数据给的是 \(U^L\)
- 那必须先做

\[
U = U^L - U^s
\]

再传给 `Uv, Uw`

否则会把 Stokes drift 在 `UL = U + Us` 中重复计算一次。

---

### 14.2 `y` 的数组方向一定要和脚本一致

当前实现中，`RectGrid` 的 `y_v/y_w/y_int` 实际是：

- 数组首点：海底 \(-H\)
- 数组末点：海表 \(0\)

如果你的上层脚本假定“首点是海表”，那么：

- 剖面插值会颠倒
- 边界条件解释会颠倒
- 表面/底部的模态与 forcing 尖峰会看起来完全错误

---

### 14.3 对 imported base profile 的导数

对于从 LES 导入的均匀网格 profile，最好不要“先插值到 Chebyshev 网格再用谱微分矩阵求导”，而应：

1. 在原始 LES 均匀网格上做有限差分得到 \(U',U'',\nu',\nu''\)
2. 再把这些导数插值到 `y_v` / `y_w`

否则很容易在 \(w\) 模式中看到沿 \(y\) 的高频伪振荡。

---

## 15. 总结

你这份实现可以概括为：

1. 用附录 B 的矩形谱配置法构造两套状态网格
2. 在 interior first-kind 节点上 enforce \(v,\omega_y\) 的 PDE
3. 在 second-kind 节点上附加边界条件
4. 构造离散输入算子 \(B\) 与输出算子 \(C\)
5. 得到离散传递矩阵

\[
T = C M^{-1} B
\]

6. 通过奇异值分解分析输入-输出放大

如果你后续要继续排查“为什么和 LES 模态对不上”，优先检查：

- 传入的是不是 \(U\) 而非 \(U^L\)
- `Us(y)` 是否与你的 LES 一致
- 垂向坐标方向是否一致
- imported profile 的导数是不是在原始均匀 LES 网格上算的
- 当前主导模态是否因底边界 no-slip 而已经偏离论文 Fig.4 的物理类型

---

## 16. 当底边界从 stress-free 改为 no-slip 时，公式和代码分别怎么变

这一节专门总结：

- 论文原始问题：上下边界均为 stress-free
- 修改后问题：顶部保持 stress-free，底部改为 no-slip

两者在连续方程、\((v,\omega_y)\) 变量形式和代码实现中分别有哪些变化。

---

### 16.1 原论文的边界条件（上下都为 stress-free）

正文 §2.1 给出的物理变量边界条件是：

\[
\frac{\partial u}{\partial y} = \frac{\partial w}{\partial y} = 0,
\qquad
v = 0
\qquad \text{at } y=0,-H.
\]

其物理含义是：

- \(v=0\)：无穿透
- \(\partial_y u = 0\)、\(\partial_y w = 0\)：无切向应力（free-slip / stress-free）

在论文使用的 \((v,\omega_y)\) 状态变量中，这等价于

\[
\hat v = 0, \qquad D^2 \hat v = 0, \qquad D \hat \omega_y = 0
\qquad \text{at } y=0,-H.
\]

也就是说：

- 对四阶的 \(v\)-equation：每个边界给两个条件
- 对二阶的 \(\omega_y\)-equation：每个边界给一个条件

总共是 \(4+2\) 个边界条件。

---

### 16.2 当底边界改为 no-slip 时，物理变量边界条件如何改变

底边界 \(y=-H\) 若改成 no-slip，则物理变量条件变为

\[
u = 0, \qquad v = 0, \qquad w = 0
\qquad \text{at } y=-H.
\]

顶部 \(y=0\) 若仍保持论文原始设定，则仍然是

\[
\frac{\partial u}{\partial y} = \frac{\partial w}{\partial y} = 0,
\qquad
v = 0
\qquad \text{at } y=0.
\]

所以整个问题变成：

- **底部**：no-slip
- **顶部**：stress-free

这是一个混合边界条件问题。

---

### 16.3 在 \((v,\omega_y)\) 变量中，底部 no-slip 等价于什么

底边界 no-slip 的三个物理条件是

\[
u=0,\qquad v=0,\qquad w=0.
\]

在 Fourier 空间中，有不可压缩条件

\[
ik_x u + Dv + ik_z w = 0.
\]

若边界上 \(u=w=0\)，则自动得到

\[
Dv = 0.
\]

此外，垂向涡量定义为

\[
\omega_y = \frac{\partial u}{\partial z} - \frac{\partial w}{\partial x}
= ik_z u - ik_x w.
\]

在边界上 \(u=w=0\) 时，就有

\[
\omega_y = 0.
\]

因此，**底边界 no-slip 在 \((v,\omega_y)\) 变量中的等价条件**是

\[
v = 0, \qquad Dv = 0, \qquad \omega_y = 0
\qquad \text{at } y=-H.
\]

这与论文原来的 stress-free 条件

\[
v = 0, \qquad D^2 v = 0, \qquad D\omega_y = 0
\qquad \text{at } y=-H
\]

相比，变化是：

- \(D^2 v = 0 \;\rightarrow\; Dv = 0\)
- \(D\omega_y = 0 \;\rightarrow\; \omega_y = 0\)

---

### 16.4 混合边界条件写成完整形式

若底部为 no-slip、顶部为 stress-free，则连续边界条件可写成：

#### 底部 \(y=-H\)

\[
v=0,\qquad Dv=0,\qquad \omega_y=0.
\]

#### 顶部 \(y=0\)

\[
v=0,\qquad D^2v=0,\qquad D\omega_y=0.
\]

这组条件仍然提供：

- \(v\)-equation 的 4 个条件
- \(\omega_y\)-equation 的 2 个条件

因此离散系统仍然是闭合的。

---

### 16.5 对应到代码里，边界条件行应该怎样改

如果使用的 second-kind 网格顺序是：

- 首点：海底 \(y=-H\)
- 末点：海表 \(y=0\)

那么边界条件在代码里应写成：

```julia
# Mixed boundary conditions:
#   bottom y = -H: no-slip      -> v = 0, Dv = 0, ω_y = 0
#   top    y =  0: stress-free -> v = 0, D²v = 0, Dω_y = 0

A[(2N + 1), 1:Nv] = Iv[1, :]          # v(-H) = 0
A[(2N + 2), 1:Nv] = Iv[end, :]        # v(0)  = 0
A[(2N + 3), 1:Nv] = g.Dv[1, :]        # Dv(-H) = 0
A[(2N + 4), 1:Nv] = g.D2v[end, :]     # D²v(0) = 0
A[(2N + 5), (Nv + 1):end] = Iw[1, :]      # ω_y(-H) = 0
A[(2N + 6), (Nv + 1):end] = g.Dw[end, :]  # Dω_y(0) = 0
```

这正对应于：

- 底部 no-slip：`v=0`, `Dv=0`, `ω_y=0`
- 顶部 stress-free：`v=0`, `D²v=0`, `Dω_y=0`

---

### 16.6 如果使用的是你这份 `build_M_B_C_rect` 风格的代码，该如何理解

在你贴出的那版 `build_M_B_C_rect` 中，边界条件通过 `M` 的若干行附加：

```julia
M[N+1, :] .= ...
M[N+2, :] .= ...
M[N+3, :] .= ...
M[N+4, :] .= ...
M[Nv+N+1, :] .= ...
M[Nv+Nw, :] .= ...
```

这里的思想和上面完全一样：

- 前 4 条给 \(v\)
- 后 2 条给 \(\omega_y\)

只不过你需要根据数组端点和物理边界的方向，确定：

- 哪个端点是海底
- 哪个端点是海表

然后把：

- 底部的 \(D^2v\) 条件改成 \(Dv\)
- 底部的 \(D\omega_y\) 条件改成 \(\omega_y\)

而顶部仍保留论文原始的 stress-free 条件。

---

### 16.7 物理影响：为什么改成 no-slip 后主模态常常会变

底部从 stress-free 改成 no-slip，不只是“换两个边界条件行”这么简单，它还会改变主导 resolvent mode 的物理类型。

常见后果包括：

- 更容易出现靠近底边界的 streak / wall mode
- 流向分量 \(u\) 更强
- 跨流和竖向分量 \(v,w\) 相对减弱
- 主导模态不再像论文 Fig.4 那样成长为贯穿全水深的对流滚涡

因此，如果你的目标是严格复现论文 Fig.4 的 full-depth Langmuir roll，
那么不应把底边界改为 no-slip。

---

### 16.8 一句话总结

底边界从 stress-free 改为 no-slip 时：

- **连续公式层面**：底部由
  \[
  \partial_y u = \partial_y w = 0,\; v=0
  \]
  改为
  \[
  u=w=v=0
  \]
- **\((v,\omega_y)\) 变量层面**：底部由
  \[
  v=0,\; D^2v=0,\; D\omega_y=0
  \]
  改为
  \[
  v=0,\; Dv=0,\; \omega_y=0
  \]
- **代码层面**：把底部边界对应的两条
  `D2v[...]` 和 `Dw[...]`
  分别替换为
  `Dv[...]` 和 `Iw[...]`

而顶部若保持论文原设定，则仍然是 stress-free，不变。

