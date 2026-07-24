# data_test

本仓库现包含 **Ocean1DRANS**：用 Julia 实现的海洋一维 RANS 水柱模型，支持 Langmuir 湍流参数化，可输出稳态背景流与湍流粘性廓线。

详见 [`Ocean1DRANS/README.md`](Ocean1DRANS/README.md)。

```bash
cd Ocean1DRANS
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. examples/xuan_shen_steady.jl
julia --project=. -e 'using Pkg; Pkg.test()'
```
