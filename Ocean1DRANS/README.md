# Ocean1DRANS

一维海洋 RANS 水柱模型（Julia）：给定外强迫 / 边界 / 初始条件，输出稳态背景流与湍流粘性，并含 Langmuir 参数化。

## 为什么初版结果对不上论文 Fig.2？

你上传的 `*_klstokes.csv` / `*_kpplt.csv` 里出现 `U~O(40)`、`ν_t~0.03`，而论文是：

| 量 | 论文 Fig.2 | 初版 klstokes |
|----|------------|---------------|
| 曲线含义 | **Lagrangian** `U_L=U+U_s` | 误把欧拉 `U` 当背景流 |
| 欧拉平均流 | 几乎可忽略 | `U/u★ ~ 40`（过大） |
| `ν_t/(u★H)` 峰值 (La=0.3) | ~0.38 | ~0.03（偏小约 10 倍） |

**根因有三：**

1. **对比错了变量**：论文 Fig.2(a) 是 `U_L`，不是欧拉 `U`。
2. **动量闭合缺 Stokes 项**：初版用 `τ=ν_t ∂U/∂z`。Langmuir 应用  
   `τ = ν_t (∂U/∂z + α_s ∂U_s/∂z)`（`α_s=1`，Lagrangian / Harcourt）。否则必须靠巨大欧拉剪切来扛应力。
3. **局部 k–ℓ 的结构性缺陷**：`k₀H=3.5` 时 Stokes 剪切只存在于近表层；局部 `P_S∝∂U_s/∂z` **无法**把中层 `ν_t` 抬到 LES 水平。要复现 Fig.2b，应使用 **`closure=:les`**（数字化 LES `ν_t`）。

## 正确用法（对接论文 / resolvent）

```julia
using Ocean1DRANS

cfg = xuan_shen_config(Nz=96, La_t=0.3, closure=:les)  # 默认已是 :les
sol = run_to_steady(cfg)

Us = sol.config.stokes.us_c
UL = sol.state.U .+ Us     # ← 对应 Fig.2(a)
νt = sol.state.νt_c        # ← 对应 Fig.2(b)

write_profiles_csv("profiles.csv", sol)  # 含 UL 列
```

```bash
julia --project=. examples/xuan_shen_steady.jl
julia --project=. examples/export_resolvent_base.jl
julia --project=. examples/validate_physics.jl
```

修正后 `closure=:les` 典型结果（La_t=0.3）：`max|U|/max|Us|≈0.13`，`max ν_t≈0.38`。

## 闭合一览

| `closure` | 说明 |
|-----------|------|
| `:les`（默认） | Fig.2b LES `ν_t` + Lagrangian 应力；论文/resolvent 推荐 |
| `:kpplt` | 峰值校准 KPPLT（`C_w≈3.6`） |
| `:klstokes` | k–ℓ + Stokes 生产；可看 Langmuir 趋势，勿直接对 Fig.2 |

动量通量一律支持 `α_s`（默认 1）。

## 物理验证

`examples/validate_physics.jl`：Lagrangian 应力平衡、Stokes/`La_t`、欧拉力≪Stokes、TKE 平衡、Langmuir 趋势、LES 峰值复现等（当前 11/11 PASS）。
