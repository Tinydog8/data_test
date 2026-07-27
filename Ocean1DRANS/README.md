# Ocean1DRANS

一维海洋 RANS 水柱模型（Julia）：给定外强迫 / 边界 / 初始条件，输出稳态背景流与湍流粘性，并含 Langmuir 参数化。

## 与成熟模型对齐

默认闭合为 **Mellor–Yamada 2.5 + Kantha–Clayson (2004)**（GOTM 同族）：

| 方程 | 形式 |
|------|------|
| ``q²`` | ``∂t q² = 2(P + P_s) - 2q³/(B1ℓ) + ∂z(Sq qℓ ∂q²/∂z)`` |
| ``q²ℓ`` | ``∂t(q²ℓ) = ℓ(E1 P + E6 P_s) - (q³/B1)(1+E2(ℓ/Lz)²) + …`` |
| ``KM`` | ``KM = q ℓ S_M``，``S_M≈0.393``，``E6=4``（KC04） |

Stokes 生产 ``P_s = KM(∂U/∂z·∂Us/∂z+…)`` 进入**两个**方程；``E6=4`` 是把 ``KM`` 抬到 LES 量级的关键（KC04 Fig.1）。严格 KC04 动量取 ``αs=0``；可选 `:harcourt` 用 Lagrangian 应力。

## 用法

```julia
using Ocean1DRANS

cfg = xuan_shen_config(Nz=96, La_t=0.3, closure=:my25)  # 默认
sol = run_to_steady(cfg)

UL = sol.state.U .+ sol.config.stokes.us_c   # Fig.2(a) 对照用 Lagrangian
νt = sol.state.νt_c
write_profiles_csv("profiles.csv", sol)
```

```bash
julia --project=. examples/xuan_shen_steady.jl
julia --project=. examples/validate_physics.jl
julia --project=. -e 'using Pkg; Pkg.test()'
```

## 闭合一览

| `closure` | 说明 |
|-----------|------|
| `:my25` / `:kc04`（默认） | MY2.5 + KC04；与 GOTM 对齐的可用物理模型 |
| `:harcourt` | 同上 + 动量 Lagrangian 应力 (`αs=1`) |
| `:les` | Fig.2b 数字化 LES `νt`；论文/resolvent 对照 |
| `:kpplt` | 峰值校准 KPPLT |
| `:klstokes` | 简化代数 k–ℓ；勿当成熟模型 |

## 与论文 Fig.2

- Fig.2(a) 画的是 **Lagrangian** ``U_L=U+U_s``，不是欧拉 ``U``。
- `:les` 精确再现数字化 LES `νt`；`:my25` 在 ``E6=4`` 下应达到同阶 ``νt``（不再差 10 倍）。
- 简化 `:klstokes` 在 ``k0H=3.5`` 时仍会低估中层混合——这是闭合层级问题，不是求解器 bug。
