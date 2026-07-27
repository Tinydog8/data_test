# Ocean1DRANS

一维海洋 RANS 水柱模型（Julia）：给定外强迫 / 边界 / 初始条件，输出稳态背景流与湍流粘性，并含 Langmuir 参数化。

## 默认：Harcourt (2015) 完整 SMC

与 GOTM `cmue_d_h15` 对齐的二阶矩闭合：

| 项 | 形式 |
|----|------|
| 动量通量 | ``⟨u'w'⟩=-qℓ(S_M ∂U/∂z + S_S ∂Us/∂z)`` |
| 涡粘 | ``KM=S_M qℓ``，``K_M^S=S_S qℓ``（独立，不是简单 ``αs``） |
| 稳定性函数 | Harcourt ARSM（含表面邻近 SPF） |
| ``q²,q²ℓ`` | KC04 型双方程，``E6=6``（H15） |

```julia
using Ocean1DRANS
cfg = xuan_shen_config(Nz=96, La_t=0.3)          # 默认 :harcourt
sol = run_to_steady(cfg)
UL = sol.state.U .+ sol.config.stokes.us_c
νt = sol.state.νt_c      # KM
νcl = sol.state.νcl_c    # K_M^S
```

## 闭合一览

| `closure` | 说明 |
|-----------|------|
| `:harcourt` / `:h15`（默认） | Harcourt 2015 完整 SMC |
| `:my25` / `:kc04` | MY2.5 + KC04（``E6=4``，``αs=0``） |
| `:les` | Fig.2b 数字化 LES `νt` |
| `:kpplt` | 峰值校准 KPPLT |
| `:klstokes` | 简化代数 k–ℓ |

对照 Fig.2 时画 **UL** 与 **nu_t**；优先看 `xuan_shen_La0.3_harcourt.csv` 与 `*_les.csv`。

## 验证 MY25 / KC04 是否实现正确

**不要用 Xuan–Shen Fig.2 判断 MY25 的 shape。** 应对照：

1. `data/kc04_fig1_KM.csv` — Kantha & Clayson (2004) Fig.1 数字化 `KM/(u★ zi)`
2. `data/mcwilliams1997_fig3b_KM.csv` — McWilliams et al. (1997) Fig.3b LES

```bash
julia --project=. examples/compare_my25_kc04.jl
# 主图：output/my25_channel_E6_compare.csv   (nu_t_E6_4 vs nu_t_E6_0)
# 叠画：output/my25_vs_kc04_fig1.csv         (列名 KM_MY25_channel_E6_*)
# 摘要：output/E6_gate_summary.txt
```

硬指标：通道应力平衡上 `E6=4` 相对 `E6=0` 必须明显抬高 `KM`（ratio > 1.5；典型 ≈2.2）。

若叠画里看起来 E6=4≈E6=0：多半画的是**旧列** `KM_MY25_E6_*`（来自未充分 spin-up 的 McWilliams+Coriolis，Ps≈0）。
请改用 `my25_channel_E6_compare.csv` 或新列名 `KM_MY25_channel_E6_*`。
