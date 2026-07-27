mutable struct ColumnState{T}
    U::Vector{T}
    V::Vector{T}
    k::Vector{T}       # TKE (= q²/2)
    ℓ::Vector{T}
    q2::Vector{T}      # twice TKE
    q2l::Vector{T}     # q²ℓ
    νt_f::Vector{T}    # KM at faces
    νt_c::Vector{T}    # KM at centers
    νcl_f::Vector{T}   # K_M^S (Stokes eddy viscosity) at faces
    νcl_c::Vector{T}   # K_M^S at centers
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
    clos = cfg.closure
    if isnothing(cfg.initial.k)
        if clos isa MY25KC04Closure || clos isa Harcourt2015Closure
            B1 = clos.B1
            k0 = 0.5 * (B1^(2 / 3)) * cfg.forcing.u★^2
        elseif clos isa KLStokesClosure
            k0 = (cfg.forcing.u★^2) / sqrt(clos.cμ)
        else
            k0 = (cfg.forcing.u★^2) / sqrt(T(0.09))
        end
        k = fill(k0, g.Nz)
    else
        k = _eval_initial(cfg.initial.k, g.zc)
    end
    κ = if clos isa MY25KC04Closure || clos isa Harcourt2015Closure
        clos.κ
    elseif clos isa KLStokesClosure
        clos.κ
    elseif clos isa KPPLTClosure
        clos.κ
    else
        0.4
    end
    ℓ_max = clos isa KLStokesClosure ? clos.ℓ_max : Inf
    channel = clos isa KLStokesClosure ? clos.channel : false
    ℓ = mixing_length(g; κ = κ, ℓ_max = ℓ_max, channel = channel)
    q2 = similar(k)
    q2l = similar(k)
    if clos isa Harcourt2015Closure
        initialize_harcourt!(q2, q2l, g, cfg.forcing.u★, clos)
        sync_tke_from_q2_h!(k, ℓ, q2, q2l, clos, g.H)
    elseif clos isa MY25KC04Closure
        initialize_my25!(q2, q2l, g, cfg.forcing.u★, clos)
        sync_tke_from_q2!(k, ℓ, q2, q2l, clos, g.H)
    else
        @. q2 = 2 * k
        @. q2l = q2 * ℓ
    end
    νt_f = zeros(T, g.Nz + 1)
    νt_c = zeros(T, g.Nz)
    νcl_f = zeros(T, g.Nz + 1)
    νcl_c = zeros(T, g.Nz)
    return ColumnState{T}(U, V, k, ℓ, q2, q2l, νt_f, νt_c, νcl_f, νcl_c, zero(T))
end

function closure_αs(clos)
    if clos isa Harcourt2015Closure
        return NaN   # not a single αs; uses independent K_M^S
    elseif hasproperty(clos, :αs)
        return clos.αs
    end
    return 0.0
end

function update_viscosity!(state::ColumnState, cfg::ModelConfig)
    clos = cfg.closure
    if clos isa Harcourt2015Closure
        # provisional SPF=1; equilibrate/advance refreshes SPF + KMS
        SPF = ones(eltype(state.q2), cfg.grid.Nz)
        Sm = similar(state.q2)
        Ss = similar(state.q2)
        Sh = similar(state.q2)
        harcourt_diffusivities!(state.νt_c, state.νcl_c, Sm, Ss, Sh,
                                state.q2, state.q2l, state.U, state.V,
                                cfg.stokes, SPF, clos, cfg.grid)
        faces_from_centers!(state.νt_f, state.νt_c, clos.νt_min)
        faces_from_centers!(state.νcl_f, state.νcl_c, clos.νt_min)
        sync_tke_from_q2_h!(state.k, state.ℓ, state.q2, state.q2l, clos, cfg.grid.H)
    elseif clos isa MY25KC04Closure
        fill!(state.νcl_c, 0)
        fill!(state.νcl_f, 0)
        eddy_viscosity_my25!(state.νt_f, state.νt_c, state.q2, state.q2l, clos, cfg.grid)
        sync_tke_from_q2!(state.k, state.ℓ, state.q2, state.q2l, clos, cfg.grid.H)
    elseif clos isa KLStokesClosure
        fill!(state.νcl_c, 0)
        fill!(state.νcl_f, 0)
        mixing_length!(state.ℓ, cfg.grid;
                       κ = clos.κ, ℓ_max = clos.ℓ_max, channel = clos.channel)
        eddy_viscosity_faces!(state.νt_f, state.k, state.ℓ, clos, cfg.grid)
        @inbounds for i in 1:cfg.grid.Nz
            state.νt_c[i] = 0.5 * (state.νt_f[i] + state.νt_f[i + 1])
        end
    elseif clos isa KPPLTClosure
        fill!(state.νcl_c, 0)
        fill!(state.νcl_f, 0)
        eddy_viscosity_kpp!(state.νt_f, clos, cfg.grid, cfg.forcing, cfg.La_t)
        @inbounds for i in 1:cfg.grid.Nz
            state.νt_c[i] = 0.5 * (state.νt_f[i] + state.νt_f[i + 1])
        end
    elseif clos isa LESNutClosure
        fill!(state.νcl_c, 0)
        fill!(state.νcl_f, 0)
        eddy_viscosity_les!(state.νt_f, clos, cfg.grid, cfg.forcing)
        @inbounds for i in 1:cfg.grid.Nz
            state.νt_c[i] = 0.5 * (state.νt_f[i] + state.νt_f[i + 1])
        end
    else
        error("unsupported closure $(typeof(clos))")
    end
    return state
end

"""
由应力平衡重建欧拉速度。

标准 / αs 形式：
```
τ = (ν+KM)(∂U/∂z + αs ∂Us/∂z)
```
完整 Harcourt：
```
τ = (ν+KM) ∂U/∂z + K_M^S ∂Us/∂z
⇒ ∂U/∂z = (τ - K_M^S ∂Us/∂z)/(ν+KM)
```
"""
function reconstruct_velocity_from_stress!(φ::AbstractVector, νt_f::AbstractVector,
                                           τ_top::Real, Fbody::Real,
                                           grid::UniformColumnGrid, ν::Real;
                                           αs::Real = 0.0,
                                           dUsdz_f = nothing,
                                           νcl_f = nothing)
    Nz = grid.Nz
    dz = grid.dz
    φ[1] = zero(eltype(φ))
    @inbounds for i in 1:Nz-1
        z_face = grid.zf[i + 1]
        νe = ν + νt_f[i + 1]
        τ = τ_top - Fbody * z_face
        if νcl_f !== nothing && dUsdz_f !== nothing
            τ_eff = τ - νcl_f[i + 1] * dUsdz_f[i + 1]
            shear = τ_eff / max(νe, eps(typeof(νe)))
        else
            shear = τ / max(νe, eps(typeof(νe)))
            if αs != 0 && !isnan(αs) && dUsdz_f !== nothing
                shear -= αs * dUsdz_f[i + 1]
            end
        end
        φ[i + 1] = φ[i] + shear * dz
    end
    return φ
end

function _reconstruct_uv!(state::ColumnState, cfg::ModelConfig)
    if cfg.closure isa Harcourt2015Closure
        reconstruct_velocity_from_stress!(state.U, state.νt_f, cfg.forcing.τx,
                                          cfg.forcing.Fx, cfg.grid, cfg.forcing.ν;
                                          dUsdz_f = cfg.stokes.dusdz_f, νcl_f = state.νcl_f)
        reconstruct_velocity_from_stress!(state.V, state.νt_f, cfg.forcing.τy,
                                          cfg.forcing.Fy, cfg.grid, cfg.forcing.ν;
                                          dUsdz_f = cfg.stokes.dvsdz_f, νcl_f = state.νcl_f)
    else
        αs = closure_αs(cfg.closure)
        reconstruct_velocity_from_stress!(state.U, state.νt_f, cfg.forcing.τx,
                                          cfg.forcing.Fx, cfg.grid, cfg.forcing.ν;
                                          αs = αs, dUsdz_f = cfg.stokes.dusdz_f)
        reconstruct_velocity_from_stress!(state.V, state.νt_f, cfg.forcing.τy,
                                          cfg.forcing.Fy, cfg.grid, cfg.forcing.ν;
                                          αs = αs, dUsdz_f = cfg.stokes.dvsdz_f)
    end
    return state
end

function equilibrate_tke!(k::AbstractVector, U::AbstractVector, V::AbstractVector,
                          stokes::StokesDrift, νt_c::AbstractVector, ℓ::AbstractVector,
                          clos::KLStokesClosure, grid::UniformColumnGrid)
    Nz = grid.Nz
    dz = grid.dz
    cε = clos.cε
    E6 = clos.E6
    αs = clos.αs
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
        P, PS = tke_production(Uz, Vz, stokes.dusdz_c[i], stokes.dvsdz_c[i],
                               νt_c[i], E6, αs)
        Prod = max(P + PS, zero(P))
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
    # Harcourt: add explicit Stokes-flux divergence ∂z(KMS ∂Us/∂z)
    if cfg.closure isa Harcourt2015Closure
        _add_stokes_flux_tendency!(state.U, state.νcl_f, cfg.stokes.dusdz_f, g, dt)
        _add_stokes_flux_tendency!(state.V, state.νcl_f, cfg.stokes.dvsdz_f, g, dt)
        Sh = fill(0.4, g.Nz)
        SPF = ones(g.Nz)
        advance_harcourt!(state.q2, state.q2l, state.U, state.V, cfg.stokes,
                          state.νt_c, state.νcl_c, Sh, SPF, cfg.closure, g, forc.u★, dt)
        update_viscosity!(state, cfg)
    elseif cfg.closure isa MY25KC04Closure
        advance_my25!(state.q2, state.q2l, state.U, state.V, cfg.stokes,
                      state.νt_c, state.νt_f, cfg.closure, g, forc.u★, dt)
        update_viscosity!(state, cfg)
    elseif cfg.closure isa KLStokesClosure
        update_tke!(state.k, state.U, state.V, cfg.stokes, state.νt_f, state.ℓ,
                    cfg.closure, g, forc, dt)
        update_viscosity!(state, cfg)
    end
    state.t += dt
    return state
end

function _add_stokes_flux_tendency!(φ::AbstractVector, νcl_f::AbstractVector,
                                    dUsdz_f::AbstractVector, grid::UniformColumnGrid, dt::Real)
    Nz = grid.Nz
    dz = grid.dz
    @inbounds for i in 1:Nz
        Ftop = νcl_f[i + 1] * dUsdz_f[i + 1]
        Fbot = νcl_f[i] * dUsdz_f[i]
        φ[i] += dt * (Ftop - Fbot) / dz
    end
    return φ
end

function estimate_dt(state::ColumnState, cfg::ModelConfig; cfl::Real = 0.3)
    νmax = maximum(cfg.forcing.ν .+ state.νt_f .+ state.νcl_f)
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
        callback !== nothing && callback(state, cfg, n)
    end
    return state
end

"""
    run_to_steady(cfg; kwargs...) -> SteadySolution

迭代求解稳态背景流与湍流粘性。无 Coriolis 时用应力平衡；
有 Coriolis 时用时间推进。
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
    αs = closure_αs(cfg.closure)
    clos_name = string(nameof(typeof(cfg.closure)))

    if verbose
        @printf("Ocean1DRANS steady run: Nz=%d, La_t=%.3f, αs=%s, closure=%s, method=%s, tol=%.1e\n",
                cfg.grid.Nz, cfg.La_t,
                isnan(αs) ? "KMS" : @sprintf("%.2f", αs),
                clos_name,
                use_stress_balance ? "stress-balance" : "time-march", tol)
    end

    if use_stress_balance
        ν_old = copy(state.νt_c)
        U_old = copy(state.U)
        k_old = copy(state.k)
        while n < max_steps
            n += 1

            if cfg.closure isa Harcourt2015Closure
                # Consistent fixed-point: reconstruct → ARSM+SPF equilibrate → diffusivities
                U_old .= state.U
                ν_old .= state.νt_c
                _reconstruct_uv!(state, cfg)
                Sh = fill(T(0.4), cfg.grid.Nz)
                equilibrate_harcourt!(state.q2, state.q2l, state.U, state.V, cfg.stokes,
                                      state.νt_c, state.νcl_c, Sh, cfg.closure, cfg.grid,
                                      cfg.forcing.u★; underrelax = min(underrelax, 0.35))
                faces_from_centers!(state.νt_f, state.νt_c, cfg.closure.νt_min)
                faces_from_centers!(state.νcl_f, state.νcl_c, cfg.closure.νt_min)
                sync_tke_from_q2_h!(state.k, state.ℓ, state.q2, state.q2l,
                                    cfg.closure, cfg.grid.H)
                _reconstruct_uv!(state, cfg)
            else
                update_viscosity!(state, cfg)
                if n > 1 && !(cfg.closure isa LESNutClosure) && !(cfg.closure isa KPPLTClosure)
                    @. state.νt_c = underrelax * state.νt_c + (1 - underrelax) * ν_old
                    state.νt_f[1] = zero(T)
                    state.νt_f[end] = zero(T)
                    @inbounds for i in 2:cfg.grid.Nz
                        state.νt_f[i] = 0.5 * (state.νt_c[i - 1] + state.νt_c[i])
                    end
                end
                _reconstruct_uv!(state, cfg)

                if cfg.closure isa MY25KC04Closure
                    equilibrate_my25!(state.q2, state.q2l, state.U, state.V, cfg.stokes,
                                      state.νt_c, cfg.closure, cfg.grid, cfg.forcing.u★;
                                      underrelax = underrelax)
                    sync_tke_from_q2!(state.k, state.ℓ, state.q2, state.q2l,
                                      cfg.closure, cfg.grid.H)
                elseif cfg.closure isa KLStokesClosure
                    @inbounds for i in 1:cfg.grid.Nz
                        state.νt_c[i] = 0.5 * (state.νt_f[i] + state.νt_f[i + 1])
                    end
                    equilibrate_tke!(state.k, state.U, state.V, cfg.stokes, state.νt_c,
                                     state.ℓ, cfg.closure, cfg.grid)
                    @. state.k = underrelax * state.k + (1 - underrelax) * k_old
                end
            end

            dU = maximum(abs, state.U .- U_old)
            dν = maximum(abs, state.νt_c .- ν_old)
            scale = max(maximum(abs, state.U), maximum(abs, cfg.stokes.us_c),
                        cfg.forcing.u★, eps(T))
            residual = max(dU, dν) / scale

            if verbose && (n % 20 == 0 || residual < tol || n == 1)
                UL = maximum(abs, state.U .+ cfg.stokes.us_c)
                @printf("  iter %5d  residual=%.3e  max|U|=%.4f  max|UL|=%.4f  max(νt)=%.4e  max(νcl)=%.4e\n",
                        n, residual, maximum(abs, state.U), UL,
                        maximum(state.νt_c), maximum(state.νcl_c))
            end

            diagnostic_done = cfg.closure isa LESNutClosure || cfg.closure isa KPPLTClosure
            if residual < tol || diagnostic_done
                if cfg.closure isa Harcourt2015Closure
                    for _ in 1:15
                        U_prev = copy(state.U)
                        _reconstruct_uv!(state, cfg)
                        Sh = fill(T(0.4), cfg.grid.Nz)
                        q2_prev = copy(state.q2)
                        equilibrate_harcourt!(state.q2, state.q2l, state.U, state.V, cfg.stokes,
                                              state.νt_c, state.νcl_c, Sh, cfg.closure, cfg.grid,
                                              cfg.forcing.u★; underrelax = 0.7)
                        faces_from_centers!(state.νt_f, state.νt_c, cfg.closure.νt_min)
                        faces_from_centers!(state.νcl_f, state.νcl_c, cfg.closure.νt_min)
                        if maximum(abs, state.q2 .- q2_prev) / max(maximum(state.q2), eps(T)) < tol &&
                           maximum(abs, state.U .- U_prev) / scale < tol
                            break
                        end
                    end
                elseif cfg.closure isa MY25KC04Closure
                    for _ in 1:20
                        update_viscosity!(state, cfg)
                        _reconstruct_uv!(state, cfg)
                        q2_prev = copy(state.q2)
                        equilibrate_my25!(state.q2, state.q2l, state.U, state.V, cfg.stokes,
                                          state.νt_c, cfg.closure, cfg.grid, cfg.forcing.u★;
                                          underrelax = 0.7)
                        if maximum(abs, state.q2 .- q2_prev) / max(maximum(state.q2), eps(T)) < tol
                            break
                        end
                    end
                elseif cfg.closure isa KLStokesClosure
                    for _ in 1:8
                        update_viscosity!(state, cfg)
                        _reconstruct_uv!(state, cfg)
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
                end
                if cfg.closure isa Harcourt2015Closure
                    faces_from_centers!(state.νt_f, state.νt_c, cfg.closure.νt_min)
                    faces_from_centers!(state.νcl_f, state.νcl_c, cfg.closure.νt_min)
                    sync_tke_from_q2_h!(state.k, state.ℓ, state.q2, state.q2l,
                                        cfg.closure, cfg.grid.H)
                else
                    update_viscosity!(state, cfg)
                end
                _reconstruct_uv!(state, cfg)
                converged = true
                residual = min(residual, tol)
                break
            end
            if !(cfg.closure isa Harcourt2015Closure)
                U_old .= state.U
                ν_old .= state.νt_c
                k_old .= state.k
            end
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
        @printf("  max|U|=%.4f  max|Us|=%.4f  max|UL|=%.4f  max(νt)=%.4e  max(νcl)=%.4e\n",
                maximum(abs, state.U), maximum(abs, cfg.stokes.us_c),
                maximum(abs, state.U .+ cfg.stokes.us_c),
                maximum(state.νt_c), maximum(state.νcl_c))
    end
    return SteadySolution(cfg, state, converged, n, T(residual), wall)
end

"""
诊断应力。Harcourt：`τ = (ν+KM)Uz + KMS Usz`；其余可用 αs 形式。
"""
function diagnostic_stress_balance(sol::SteadySolution)
    cfg = sol.config
    g = cfg.grid
    νe = cfg.forcing.ν .+ sol.state.νt_f
    αs = closure_αs(cfg.closure)
    τ_f = similar(g.zf)
    Uz_f = similar(g.zf)
    @inbounds for i in eachindex(g.zf)
        τ_f[i] = cfg.forcing.τx - cfg.forcing.Fx * g.zf[i]
        if cfg.closure isa Harcourt2015Closure
            Uz_f[i] = (τ_f[i] - sol.state.νcl_f[i] * cfg.stokes.dusdz_f[i]) /
                      max(νe[i], eps(eltype(νe)))
        else
            a = isnan(αs) ? 0.0 : αs
            Uz_f[i] = τ_f[i] / max(νe[i], eps(eltype(νe))) - a * cfg.stokes.dusdz_f[i]
        end
    end
    return (; νe, τ_f, Uz_f, αs, νcl = sol.state.νcl_f)
end
