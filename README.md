# data_test

本仓库包含 **Ocean1DRANS**：Julia 一维海洋 RANS（含 Langmuir），用于输出稳态背景流与湍流粘性。

**与 Xuan & Shen (2025) 对照时**：请用 `closure=:les`，并比较 CSV 中的 **`UL`**（不是欧拉 `U`）。详见 [`Ocean1DRANS/README.md`](Ocean1DRANS/README.md)。

```bash
cd Ocean1DRANS
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. examples/xuan_shen_steady.jl
julia --project=. examples/validate_physics.jl
julia --project=. -e 'using Pkg; Pkg.test()'
```
