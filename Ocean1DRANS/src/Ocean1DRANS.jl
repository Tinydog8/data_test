"""
    Ocean1DRANS

一维海洋 RANS 水柱模型，支持 Langmuir 湍流参数化，
用于在给定外强迫、边界条件与初始条件下求解稳态背景流与湍流粘性廓线。

物理参考：
- Harcourt (2015) 完整 SMC（GOTM `cmue_d_h15`）：ARSM 稳定性函数 ``S_M,S_S`` + SPF
- GOTM / Mellor–Yamada 2.5 + Kantha & Clayson (2004): 预后 q²、q²ℓ，E6 Stokes 源
- McWilliams et al. (1997) / Li & Fox-Kemper (2017): KPPLT
- Xuan & Shen (2025) 型无分层 Langmuir 通道设定（压力梯度平衡风应力）
"""
module Ocean1DRANS

using LinearAlgebra
using Printf

include("utils.jl")
include("grid.jl")
include("stokes.jl")
include("config.jl")
include("closures/mixing_length.jl")
include("closures/kl_stokes.jl")
include("closures/my25_kc04.jl")
include("closures/harcourt_smc.jl")
include("closures/kpp_lt.jl")
include("closures/les_nut.jl")
include("presets.jl")
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
    Harcourt2015Closure, HarcourtMomentumClosure,
    MY25KC04Closure, KanthaClayson2004Closure, KC04LagrangianClosure,
    KLStokesClosure, KPPLTClosure, LESNutClosure,
    # solve
    ColumnState, SteadySolution,
    integrate!, run_to_steady,
    diagnostic_stress_balance,
    # io / utils
    write_profiles_csv, profile_dict,
    load_les_nut_csv,
    # validation
    CheckResult, run_physics_validation

end # module
