"""
    Ocean1DRANS

一维海洋 RANS 水柱模型，支持 Langmuir 湍流参数化，
用于在给定外强迫、边界条件与初始条件下求解稳态背景流与湍流粘性廓线。

物理参考：
- GOTM / Kantha & Clayson (2004): TKE 中的 Stokes 剪切生产
- Harcourt (2013, 2015): Langmuir second-moment closure 思想
- McWilliams et al. (1997) / Li & Fox-Kemper (2017): KPPLT
- Xuan & Shen (2025) 型无分层 Langmuir 通道设定（压力梯度平衡风应力）
"""
module Ocean1DRANS

using LinearAlgebra
using Printf

include("grid.jl")
include("stokes.jl")
include("config.jl")
include("closures/mixing_length.jl")
include("closures/kl_stokes.jl")
include("closures/kpp_lt.jl")
include("momentum.jl")
include("solver.jl")
include("io.jl")
include("validation.jl")

export
    # grid
    UniformColumnGrid,
    cell_centers, interfaces, depths,
    # Stokes
    StokesDrift,
    monochromatic_stokes, exponential_stokes,
    # config
    Forcing, BoundarySetup, InitialState, ModelConfig,
    xuan_shen_config, mcwilliams1997_config,
    # closures
    KLStokesClosure, KPPLTClosure,
    # solve
    ColumnState, SteadySolution,
    integrate!, run_to_steady,
    diagnostic_stress_balance,
    # io
    write_profiles_csv, profile_dict,
    # validation
    CheckResult, run_physics_validation

end # module
