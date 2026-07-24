# Ocean1DRANS

一维海洋 RANS 水柱模型（Julia），在给定外强迫、边界条件与初始条件下求解**稳态背景流**与**湍流粘性系数廓线**，并包含 **Langmuir 湍流参数化**。

## 物理框架

动量方程（水平均匀水柱）：

```text
∂U/∂t =  f (V + Vˢ) + Fₓ + ∂z[(ν + νₜ) ∂z U]
∂V/∂t = -f (U + Uˢ) + Fᵧ + ∂z[(ν + νₜ) ∂z V]
```

其中 `(Uˢ, Vˢ)` 为 Stokes 漂移，`f(Vˢ, -Uˢ)` 为 Stokes–Coriolis 力。

### Langmuir 参数化（两种闭合）

1. **`KLStokesClosure`**（默认）— 参考 GOTM / Kantha & Clayson (2004)
   - 预后 TKE，代数混合长度
   - TKE 源项含欧拉剪切生产 `P` 与 Stokes 剪切生产 `E₆ P_S`
   - `νₜ = cμ √k · ℓ`

2. **`KPPLTClosure`** — 参考 Large et al. (1994) + Li & Fox-Kemper (2017)
   - K 廓线形状 `G(σ)=σ(1-σ)²`
   - Langmuir 增强速度尺度 `wₛ ← wₛ √(1 + C_w / La_t²)`

Stokes 漂移默认深水指数型（McWilliams et al. 1997；Xuan & Shen 2025）：

```text
Uˢ(z) = Uˢ₀ exp(2 k₀ z),   La_t = √(u★ / Uˢ₀),   Uˢ₀ = u★ / La_t²
```

## 快速开始

```julia
using Ocean1DRANS

# Xuan & Shen 型无分层 Langmuir 通道（无量纲：u★=1, H=1, Reτ=1000）
cfg = xuan_shen_config(Nz=96, La_t=0.3, Reτ=1000, k0H=3.5, E6=4.0)
sol = run_to_steady(cfg; tol=2e-5, verbose=true)

# 稳态廓线
z   = sol.config.grid.zc
U   = sol.state.U
νt  = sol.state.νt_c
Us  = sol.config.stokes.us_c

write_profiles_csv("profiles.csv", sol)
```

预设开洋混合层（含 Coriolis）：

```julia
cfg = mcwilliams1997_config(Nz=64, H=90.0, u★=0.0061, La_t=0.3)
sol = run_to_steady(cfg)
```

## 运行算例 / 测试

```bash
cd Ocean1DRANS
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. examples/xuan_shen_steady.jl
julia --project=. examples/mcwilliams1997_steady.jl
julia --project=. examples/validate_physics.jl   # 物理准确性验证（含 LES 对照）
julia --project=. -e 'using Pkg; Pkg.test()'
```

### 物理验证内容（`validate_physics.jl`）

| 检查 | 含义 |
|------|------|
| 应力平衡 | 底应力 0、表应力 `u★²`、离散 `νe ∂U/∂z` 重建 |
| Stokes / `La_t` | 指数廓线与 `La_t=√(u★/Us0)` |
| 剪切结构 | 风生通道 `∂U/∂z≥0` |
| TKE 平衡 | `P + E₆ P_S ≈ ε`（k–ℓ 局部平衡） |
| Langmuir 趋势 | `E₆`、更小 `La_t`、KPPLT 开关均增强 `νₜ` |
| KPPLT 峰 | `G(σ)` 峰值位于 `σ≈1/3` |
| LES 形态 | 对照 `data/les_eddy_viscosity_fig2b.csv`（Xuan & Shen Fig.2b）的近壁衰减与内部单峰 |

## 与本仓库其它工作的衔接

本模块输出的 `(U(z), νₜ(z), Uˢ(z))` 可直接作为 Langmuir resolvent / 线性稳定性分析的基流与涡粘剖面输入，替代此前从 LES 数字化的 `nu_t` 剖面（参见仓库中 `langmuir_resolvent_cpu.jl` 的 `nuT_profile` 接口）。

```bash
julia --project=. examples/export_resolvent_base.jl
```

生成 `output/resolvent_base_La*_*.csv`，列为 `y/H, U/u★, Us/u★, νt/(u★H)`。

可通过增大 `KPPLTClosure(Cw=...)` 或 `KLStokesClosure(E6=...)` 校准涡粘量级，使其接近 LES（例如 Xuan & Shen Fig.2b 中 `La_t=0.3` 时 `νt/(u★H)` 峰值约 0.38）。

## 主要参考

- Kantha & Clayson (2004), *On the effect of surface gravity waves on mixing in an oceanic mixed layer*
- Harcourt (2013, 2015), Langmuir second-moment closures
- Li & Fox-Kemper (2017), KPPLT entrainment enhancement
- McWilliams et al. (1997), Langmuir turbulence LES
- GOTM / CVMix Langmuir modules
- OceanTurb.jl（Julia 一维海洋湍流参数化框架）
