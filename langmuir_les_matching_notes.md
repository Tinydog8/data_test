## 使 Oceananigans 脚本尽量匹配 Xuan & Shen (2025) LES 的修改要点

这篇论文的 LES 物理设定，与当前脚本的差异非常大。若目标是尽量复现实验设定，而不是保留现有的分层/热盐/底拖曳场景，建议按下面的原则改：

### 1. 去掉温盐与浮力

论文第 2.2 节使用的是无分层、无热盐的 canonical Langmuir turbulence。

- 删除 `tracers = (:T, :S)`
- 删除 `buoyancy = SeawaterBuoyancy()`
- 删除 `T_bcs`、`S_bcs`、`Qᵀ`、`dTdz` 等热盐相关设置
- 删除 `Tᵢ`、`Sᵢ` 初始条件

也就是说，这里应只求解速度场，不求解温盐。

### 2. 域尺寸要改成论文的长宽比

论文采用

- 论文坐标中：`Lx = 8πH`、竖直方向 `Ly = H`、跨风向 `Lz = 4πH`
- 在 Oceananigans 里通常把竖直方向放在 `z`，因此脚本中映射为：
  - `x = 8πH`
  - `y = 4πH`
  - `z = H`

若你继续令 `H = 25 m`，则应改成

- `Lx = 8π * 25 ≈ 628.3 m`
- `Ly = 4π * 25 ≈ 314.2 m`
- `Lz = 25 m`

当前脚本 `extent = (600, 600, H)` 的横向长宽比不对。

### 3. Stokes 漂移必须改成论文使用的深水指数型

论文第 2.2 节给出

`Uˢ(z) = Uˢ₀ exp(2 k₀ z)`, `k₀ H = 3.5`

因此：

- 不要再用当前代码里的有限水深双曲函数形式
- 应改成深水单色波指数衰减形式

对应到 Oceananigans 中可写成

```julia
k₀ = 3.5 / H
uˢ(z) = Uˢ₀ * exp(2k₀ * z)
∂z_uˢ(z, t) = 2k₀ * Uˢ₀ * exp(2k₀ * z)
```

### 4. 参数应由 `La_t` 与 `Reτ` 决定，而不是由经验风速和波高直接决定

论文使用两组参数：

- `La_t = 0.2` 或 `0.3`
- `Reτ = u★ H / ν = 1000`

因此更合理的做法是先选：

```julia
La_t = 0.2
u★ = 0.01
ν = u★ * H / 1000
Uˢ₀ = u★ / La_t^2
```

注意：这里的 `ν` 是为了满足论文无量纲控制参数而选的有效分子黏性，不是海水真实分子黏性。

### 5. 顶部风应力保留，但底部拖曳必须去掉

论文边界条件的关键点是：

- 表面施加常风应力
- 底部 stress-free
- 不使用对数底边界层阻力

因此：

- 删除 `drag_u`、`drag_v`
- 删除底部 roughness length、`cᴰ`、`drag_bc_u`、`drag_bc_v`
- 改成底部零切应力

Oceananigans 近似可写成

```julia
τw = -u★^2
u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(τw),
                                bottom = FluxBoundaryCondition(0.0))
v_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0.0),
                                bottom = FluxBoundaryCondition(0.0))
```

### 6. 压力梯度需要和表面风应力平衡

论文指出使用 constant adverse pressure gradient 来平衡整体动量，因此比你现在人为指定 `up = 0.01` 更合适的写法是

```julia
Fᵤ(x, y, z, t) = -u★^2 / H
```

也就是用 `u★` 决定体力项，而不是另设一个独立的 `up`。

### 7. SGS 模型应改为动态 Smagorinsky

论文明确使用 dynamic Smagorinsky model。

因此：

- 不要使用 `AnisotropicMinimumDissipation()`
- 尽量改成 `DynamicSmagorinsky()`
- 若要保留分子黏性，可与 `ScalarDiffusivity(ν=ν)` 组合

Oceananigans 中建议写成

```julia
closure = (ScalarDiffusivity(ν=ν), DynamicSmagorinsky())
```

> 注意：`DynamicSmagorinsky()` 在较新版本 Oceananigans 中才比较稳定，建议使用包含相关修复的版本。

### 8. 垂向网格要在表面和底部加密

论文使用：

- `Nx × Ny × Nz = 768 × 256 × 1024`（其中竖直方向为 256）
- 表面与底部附近网格聚簇
- `Δy⁺_min = 0.25`

Oceananigans 虽不能完全复制论文的 hybrid pseudo-spectral / finite-difference 离散，但至少应：

- 用非均匀垂向网格
- 在表面和底部加密
- 水平方向保持周期边界

### 9. 尽量只保留论文真正需要的输出

若目标是和论文第 2.3 节及后续 resolvent 分析对接，应优先输出：

- `U_L(z) = <u_L>_{x,y}`
- `U_E(z) = U_L - Uˢ`
- `u'w'(z)`（或按你坐标定义对应的垂向 Reynolds shear stress）
- `νₑ(z)` 或可用于构造 `ν_t(z)` 的量

而你当前脚本里大量高阶张量、热通量、盐度、应变率与底边界层诊断，不是复现本文 LES 的必要项。

### 10. 需要接受的一个事实：Oceananigans 只能“物理设定接近”，不能“数值格式完全相同”

论文原始 LES 使用：

- hybrid pseudo-spectral / finite-difference
- fractional-step
- dynamic Smagorinsky

而 Oceananigans 的离散格式不同，所以最合理的目标是：

1. 先把**物理设定**对齐；
2. 再用统计量（均速、Reynolds stress、Stokes 剪切、主导尺度）去比较；
3. 不要期待逐点场完全一致。

---

仓库里新增的 `langmuir_les_xuan2025_oceananigans.jl` 就是按上述原则整理的一份 Oceananigans 版本模板。
