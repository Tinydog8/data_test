mutable struct ColumnState{T}
    U::Vector{T}
    V::Vector{T}
    k::Vector{T}
    ℓ::Vector{T}
    νt_f::Vector{T}   # faces, length Nz+1
    νt_c::Vector{T}   # centers
    t::T
end

struct SteadySolution{T,C}
    config::ModelConfig{T,C}
    state::ColumnState{T}
    converged::Bool
    iterations::Int
    residual::T
    walltime::Float64
end

function _eval_initial(val, z::AbstractVector{T}) where {T}
    if val isa Number
        return fill(T(val), length(z))
    elseif val isa AbstractVector
        length(val) == length(z) || throw(ArgumentError("IC vector length mismatch"))
        return T.(collect(val))
    else
        return T[val(zi) for zi in z]
    end
end

function init_state(cfg::ModelConfig{T}) where {T}
    g = cfg.grid
    U = _eval_initial(cfg.initial.U, g.zc)
    V = _eval_initial(cfg.initial.V, g.zc)
    if isnothing(cfg.initial.k)
        cμ = cfg.closure isa KLStokesClosure ? cfg.closure.cμ : T(0.09)
        k0 = (cfg.forcing.u★^2) / sqrt(cμ)
        k = fill(k0, g.Nz)
    else
        k = _eval_initial(cfg.initial.k, g.zc)
    end
    ℓ = mixing_length(g;
                      κ = cfg.closure isa KLStokesClosure ? cfg.closure.κ :
                          cfg.closure isa KPPLTClosure ? cfg.closure.κ : 0.4,
                      ℓ_max = cfg.closure isa KLStokesClosure ? cfg.closure.ℓ_max : Inf,
                      channel = cfg.closure isa KLStokesClosure ? cfg.closure.channel : false)
    νt_f = zeros(T, g.Nz + 1)
    νt_c = zeros(T, g.Nz)
    return ColumnState{T}(U, V, k, ℓ, νt_f, νt_c, zero(T))
end

function update_viscosity!(state::ColumnState, cfg::ModelConfig)
    clos = cfg.closure
    if clos isa KLStokesClosure
        mixing_length!(state.ℓ, cfg.grid;
                       κ = clos.κ, ℓ_max = clos.ℓ_max, channel = clos.channel)
        eddy_viscosity_faces!(state.νt_f, state.k, state.ℓ, clos, cfg.grid)
    elseif clos isa KPPLTClosure
        eddy_viscosity_kpp!(state.νt_f, clos, cfg.grid, cfg.forcing, cfg.La_t)
    else
        error("unsupported closure $(typeof(clos))")
    end
    @inbounds for i in 1:cfg.grid.Nz
        state.νt_c[i] = 0.5 * (state.νt_f[i] + state.νt_f[i + 1])
    end
    return state
end

"""
无 Coriolis 时由应力平衡重建速度：
`(ν+νt) ∂φ/∂z = τ_top - F * z`（界面），再从底到表积分。
底部 stress-free 要求 `τ_top + F*H = 0`。
速度零点取底部第一层 `φ[1]=0`（规范条件）。
"""
function reconstruct_velocity_from_stress!(φ::AbstractVector, νt_f::AbstractVector,
                                           τ_top::Real, Fbody::Real,
                                           grid::UniformColumnGrid, ν::Real)
    Nz = grid.Nz
    dz = grid.dz
    φ[1] = zero(eltype(φ))
    # face shear stress and integrate cell-to-cell
    # Uz at face i+1/2 between cells i and i+1:
    @inbounds for i in 1:Nz-1
        z_face = grid.zf[i + 1]          # interface between cell i and i+1
        νe = ν + νt_f[i + 1]
        τ = τ_top - Fbody * z_face
        dφ = τ / max(νe, eps(typeof(νe))) * dz
        φ[i + 1] = φ[i] + dφ
    end
    return φ
end

"""
局部平衡 TKE：`P + E6 P_S = ε = cε k^{3/2}/ℓ`，忽略输运。
"""
function equilibrate_tke!(k::AbstractVector, U::AbstractVector, V::AbstractVector,
                          stokes::StokesDrift, νt_c::AbstractVector, ℓ::AbstractVector,
                          clos::KLStokesClosure, grid::UniformColumnGrid)
    Nz = grid.Nz
    dz = grid.dz
    cε = clos.cε
    E6 = clos.E6
    k_min = clos.k_min
    @inbounds for i in 1:Nz
        if i == 1
            Uz = (U[2] - U[1]) / dz
            Vz = (V[2] - V[1]) / dz
        elseif i == Nz
            Uz = (U[Nz] - U[Nz - 1]) / dz
            Vz = (V[Nz] - V[Nz - 1]) / dz
        else
            Uz = (U[i + 1] - U[i - 1]) / (2dz)
            Vz = (V[i + 1] - V[i - 1]) / (2dz)
        end
        P, PS = tke_production(Uz, Vz, stokes.dusdz_c[i], stokes.dvsdz_c[i], νt_c[i], E6)
        Prod = max(P + PS, zero(P))
        # ε = cε k^{3/2}/ℓ = Prod  ⇒  k = (Prod * ℓ / cε)^{2/3}
        k[i] = max((Prod * ℓ[i] / cε)^(2 / 3), k_min)
    end
    return k
end

function timestep!(state::ColumnState, cfg::ModelConfig, dt::Real)
    g = cfg.grid
    forc = cfg.forcing
    update_viscosity!(state, cfg)
    apply_body_forces!(state.U, state.V, cfg.stokes, forc, dt)
    νe = _face_viscosity(state.νt_f, forc.ν)
    diffuse_velocity!(state.U, νe, forc.τx, g, cfg.boundary, dt;
                      φ_bottom_companion = state.V)
    diffuse_velocity!(state.V, νe, forc.τy, g, cfg.boundary, dt;
                      φ_bottom_companion = state.U)
    if cfg.closure isa KLStokesClosure
        update_tke!(state.k, state.U, state.V, cfg.stokes, state.νt_f, state.ℓ,
                    cfg.closure, g, forc, dt)
        update_viscosity!(state, cfg)
    end
    state.t += dt
    return state
end

function estimate_dt(state::ColumnState, cfg::ModelConfig; cfl::Real = 0.3)
    νmax = maximum(cfg.forcing.ν .+ state.νt_f)
    dt_diff = cfl * cfg.grid.dz^2 / max(νmax, eps(typeof(νmax)))
    dt_fric = 0.1 * cfg.grid.H / max(cfg.forcing.u★, eps(typeof(cfg.forcing.u★)))
    return min(dt_diff, dt_fric)
end

function integrate!(state::ColumnState, cfg::ModelConfig;
                    dt = nothing, nsteps::Integer = 1000, callback = nothing)
    update_viscosity!(state, cfg)
    Δt = isnothing(dt) ? estimate_dt(state, cfg) : dt
    for n in 1:nsteps
        timestep!(state, cfg, Δt)
        if callback !== nothing
            callback(state, cfg, n)
        end
    end
    return state
end

"""
    run_to_steady(cfg; kwargs...) -> SteadySolution

迭代求解稳态背景流与湍流粘性廓线。

对无 Coriolis 情形（如 Xuan–Shen 通道）采用：
1. 由闭合更新 `νt`
2. 应力平衡重建 `(U,V)`
3. 对 `KLStokesClosure` 用局部生产–耗散平衡更新 `k`

对有 Coriolis 情形回退到隐式扩散时间推进。
"""
function run_to_steady(cfg::ModelConfig{T};
                       tol::Real = 1e-6,
                       max_steps::Integer = 5_000,
                       dt = nothing,
                       check_every::Integer = 50,
                       cfl::Real = 0.25,
                       verbose::Bool = true,
                       underrelax::Real = 0.5) where {T}
    t0 = time()
    state = init_state(cfg)
    update_viscosity!(state, cfg)

    residual = typemax(T)
    converged = false
    n = 0
    use_stress_balance = abs(cfg.forcing.f) < eps(T)

    if verbose
        @printf("Ocean1DRANS steady run: Nz=%d, La_t=%.3f, method=%s, tol=%.1e\n",
                cfg.grid.Nz, cfg.La_t,
                use_stress_balance ? "stress-balance" : "time-march", tol)
    end

    if use_stress_balance
        ν_old = copy(state.νt_c)
        U_old = copy(state.U)
        k_old = copy(state.k)
        while n < max_steps
            n += 1
            update_viscosity!(state, cfg)

            if n > 1
                @. state.νt_c = underrelax * state.νt_c + (1 - underrelax) * ν_old
                state.νt_f[1] = zero(T)
                state.νt_f[end] = zero(T)
                @inbounds for i in 2:cfg.grid.Nz
                    state.νt_f[i] = 0.5 * (state.νt_c[i - 1] + state.νt_c[i])
                end
            end

            reconstruct_velocity_from_stress!(state.U, state.νt_f, cfg.forcing.τx,
                                              cfg.forcing.Fx, cfg.grid, cfg.forcing.ν)
            reconstruct_velocity_from_stress!(state.V, state.νt_f, cfg.forcing.τy,
                                              cfg.forcing.Fy, cfg.grid, cfg.forcing.ν)

            if cfg.closure isa KLStokesClosure
                @inbounds for i in 1:cfg.grid.Nz
                    state.νt_c[i] = 0.5 * (state.νt_f[i] + state.νt_f[i + 1])
                end
                equilibrate_tke!(state.k, state.U, state.V, cfg.stokes, state.νt_c,
                                 state.ℓ, cfg.closure, cfg.grid)
                @. state.k = underrelax * state.k + (1 - underrelax) * k_old
            end

            dU = maximum(abs, state.U .- U_old)
            dν = maximum(abs, state.νt_c .- ν_old)
            scale = max(maximum(abs, state.U), cfg.forcing.u★, eps(T))
            residual = max(dU, dν) / scale

            if verbose && (n % 20 == 0 || residual < tol || n == 1)
                @printf("  iter %5d  residual=%.3e  max|U|=%.4f  max(νt)=%.4e\n",
                        n, residual, maximum(abs, state.U), maximum(state.νt_c))
            end

            if residual < tol
                # 最终一致性：νt → U → k → νt
                if cfg.closure isa KLStokesClosure
                    for _ in 1:5
                        update_viscosity!(state, cfg)
                        reconstruct_velocity_from_stress!(state.U, state.νt_f, cfg.forcing.τx,
                                                          cfg.forcing.Fx, cfg.grid, cfg.forcing.ν)
                        reconstruct_velocity_from_stress!(state.V, state.νt_f, cfg.forcing.τy,
                                                          cfg.forcing.Fy, cfg.grid, cfg.forcing.ν)
                        @inbounds for i in 1:cfg.grid.Nz
                            state.νt_c[i] = 0.5 * (state.νt_f[i] + state.νt_f[i + 1])
                        end
                        k_prev = copy(state.k)
                        equilibrate_tke!(state.k, state.U, state.V, cfg.stokes, state.νt_c,
                                         state.ℓ, cfg.closure, cfg.grid)
                        if maximum(abs, state.k .- k_prev) / max(maximum(state.k), eps(T)) < tol
                            break
                        end
                    end
                    update_viscosity!(state, cfg)
                    reconstruct_velocity_from_stress!(state.U, state.νt_f, cfg.forcing.τx,
                                                      cfg.forcing.Fx, cfg.grid, cfg.forcing.ν)
                    reconstruct_velocity_from_stress!(state.V, state.νt_f, cfg.forcing.τy,
                                                      cfg.forcing.Fy, cfg.grid, cfg.forcing.ν)
                else
                    update_viscosity!(state, cfg)
                    reconstruct_velocity_from_stress!(state.U, state.νt_f, cfg.forcing.τx,
                                                      cfg.forcing.Fx, cfg.grid, cfg.forcing.ν)
                    reconstruct_velocity_from_stress!(state.V, state.νt_f, cfg.forcing.τy,
                                                      cfg.forcing.Fy, cfg.grid, cfg.forcing.ν)
                end
                converged = true
                break
            end
            U_old .= state.U
            ν_old .= state.νt_c
            k_old .= state.k
        end
    else
        Δt = isnothing(dt) ? estimate_dt(state, cfg; cfl = cfl) : T(dt)
        U_old = copy(state.U)
        V_old = copy(state.V)
        ν_old = copy(state.νt_c)
        while n < max_steps
            for _ in 1:check_every
                timestep!(state, cfg, Δt)
                n += 1
                n >= max_steps && break
            end
            dU = maximum(abs, state.U .- U_old)
            dV = maximum(abs, state.V .- V_old)
            dν = maximum(abs, state.νt_c .- ν_old)
            scale = max(maximum(abs, state.U), maximum(abs, state.V),
                        cfg.forcing.u★, eps(T))
            residual = max(dU, dV, dν) / scale
            if verbose && (n % (check_every * 5) == 0 || residual < tol)
                @printf("  step %7d  t=%.3e  residual=%.3e  max|U|=%.4f  max(νt)=%.4e\n",
                        n, state.t, residual, maximum(abs, state.U), maximum(state.νt_c))
            end
            if residual < tol
                converged = true
                break
            end
            U_old .= state.U
            V_old .= state.V
            ν_old .= state.νt_c
        end
    end

    wall = time() - t0
    if verbose
        @printf("Finished: converged=%s, iters=%d, residual=%.3e, wall=%.2fs\n",
                converged, n, residual, wall)
    end
    return SteadySolution(cfg, state, converged, n, T(residual), wall)
end

"""
由稳态应力平衡诊断界面剪切：`νe ∂U/∂z = τx - Fx*z`。
"""
function diagnostic_stress_balance(sol::SteadySolution)
    cfg = sol.config
    g = cfg.grid
    νe = cfg.forcing.ν .+ sol.state.νt_f
    τ_f = similar(g.zf)
    @inbounds for i in eachindex(g.zf)
        τ_f[i] = cfg.forcing.τx - cfg.forcing.Fx * g.zf[i]
    end
    return (; νe, τ_f, Uz_f = τ_f ./ max.(νe, eps(eltype(νe))))
end
