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
| `:my25` / `:kc04` | MY2.5 + KC04（通道默认 ``E6=4``；McWilliams 默认 ``E6=7.2``） |
| `:les` | Fig.2b 数字化 LES `νt` |
| `:kpplt` | 峰值校准 KPPLT |
| `:klstokes` | 简化代数 k–ℓ |

对照 Fig.2 时画 **UL** 与 **nu_t**；优先看 `xuan_shen_La0.3_harcourt.csv` 与 `*_les.csv`。

## 验证 MY25 / KC04 Fig.1

```bash
julia --project=. examples/compare_my25_kc04.jl
# 叠画：output/my25_vs_kc04_fig1.csv
```

McWilliams / KC04 设定要点：

- `zi = H = 33 m`，稳态 **Ekman–Stokes**（不是未收敛的时间推进）
- 仅表面壁面律 + 预后 ``q²/q²ℓ``（代数局部平衡会把 ``ℓ`` 压死）
- **E6=7.2**（Kantha et al. 2010 对原文 E6=4 笔误的更正）对齐 Fig.1 粗红线量级
- `KM_MY25_noLC` / `KM_MY25_E6_4` / `KM_MY25_E6_7p2` 与数字化 `KM_KC04_*` 同轴比较

通道算例仍可用 `E6=4` 做实现硬指标（`my25_channel_E6_compare.csv`）。
