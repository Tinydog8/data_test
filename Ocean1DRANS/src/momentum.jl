"""
一维动量方程的隐式扩散 + 显式 Coriolis / Stokes–Coriolis / 体力。

```
∂U/∂t = + f (V + Vs) + Fx + ∂z[(ν+νt) ∂z U]
∂V/∂t = - f (U + Us) + Fy + ∂z[(ν+νt) ∂z V]
```

表面通量边界：`(ν+νt) ∂z U|_s = τx`（`τx` 通常为 `-u★²`）。
底部：`:stress_free` / `:no_slip` / `:quadratic_drag`。
"""

function _face_viscosity(νt_f::AbstractVector, ν::Real)
    return ν .+ νt_f
end

"""
组装隐式三对角扩散算子，施加表面通量与底部 BC，前进一步。
对 U、V 分别调用。
"""
function diffuse_velocity!(φ::AbstractVector, νe_f::AbstractVector, τ_top::Real,
                           grid::UniformColumnGrid, boundary::BoundarySetup,
                           dt::Real; φ_bottom_companion = nothing)
    Nz = grid.Nz
    dz = grid.dz
    a = zeros(eltype(φ), Nz)  # lower
    b = zeros(eltype(φ), Nz)  # diag
    c = zeros(eltype(φ), Nz)  # upper
    rhs = copy(φ)

    @inbounds for i in 1:Nz
        Km = νe_f[i]
        Kp = νe_f[i+1]
        a[i] = -dt * Km / dz^2
        c[i] = -dt * Kp / dz^2
        b[i] = 1 - a[i] - c[i]
    end

    # Top face i = Nz+1: flux condition K ∂φ/∂z = τ_top
    # Discretize: Kp (φ_ghost - φ_Nz)/dz = τ_top → eliminate ghost
    # Contribution of top flux to cell Nz: (τ_top - Km(φ_Nz-φ_{Nz-1})/dz)/dz
    # With implicit: replace c[Nz]=0 and add τ_top*dt/dz to rhs; adjust b
    c[Nz] = 0
    b[Nz] = 1 + dt * νe_f[Nz] / dz^2   # only bottom face of top cell
    rhs[Nz] += dt * τ_top / dz

    # Bottom BC
    if boundary.bottom === :stress_free
        # K ∂φ/∂z = 0 at bottom → a[1]=0
        a[1] = 0
        b[1] = 1 + dt * νe_f[2] / dz^2
    elseif boundary.bottom === :no_slip
        # φ_bottom_face = 0 → ghost such that (φ_1+φ_g)/2=0 → φ_g=-φ_1
        # flux at bottom ≈ νe (φ_1 - φ_g)/dz = 2 νe φ_1/dz
        a[1] = 0
        b[1] = 1 + dt * (νe_f[2] / dz^2 + 2νe_f[1] / dz^2)
    elseif boundary.bottom === :quadratic_drag
        # τ_b = -Cd |u| u  (explicitized with companion speed)
        speed = if isnothing(φ_bottom_companion)
            abs(φ[1])
        else
            sqrt(φ[1]^2 + φ_bottom_companion[1]^2)
        end
        Cd = oftype(dz, boundary.Cd)
        # K ∂φ/∂z|_b = τ_b = -Cd |u| φ
        # implicit linearization: τ_b ≈ -Cd*speed*φ_new
        a[1] = 0
        drag = Cd * max(speed, eps(typeof(speed)))
        b[1] = 1 + dt * νe_f[2] / dz^2 + dt * drag / dz
    else
        error("unknown bottom BC: $(boundary.bottom)")
    end

    thomas_solve!(φ, a, b, c, rhs)
    return φ
end

function thomas_solve!(x::AbstractVector, a::AbstractVector, b::AbstractVector,
                       c::AbstractVector, d::AbstractVector)
    n = length(b)
    cp = similar(c)
    dp = similar(d)
    cp[1] = c[1] / b[1]
    dp[1] = d[1] / b[1]
    @inbounds for i in 2:n
        denom = b[i] - a[i] * cp[i-1]
        cp[i] = i == n ? zero(eltype(c)) : c[i] / denom
        dp[i] = (d[i] - a[i] * dp[i-1]) / denom
    end
    x[n] = dp[n]
    @inbounds for i in n-1:-1:1
        x[i] = dp[i] - cp[i] * x[i+1]
    end
    return x
end

"""
显式 Coriolis + Stokes–Coriolis + 常体力。
"""
function apply_body_forces!(U::AbstractVector, V::AbstractVector,
                            stokes::StokesDrift, forcing::Forcing, dt::Real)
    f = forcing.f
    Fx = forcing.Fx
    Fy = forcing.Fy
    if f == 0 && Fx == 0 && Fy == 0
        return U, V
    end
    @inbounds for i in eachindex(U)
        # ∂U/∂t = f(V + Vs) + Fx
        # ∂V/∂t = -f(U + Us) + Fy
        Us = stokes.us_c[i]
        Vs = stokes.vs_c[i]
        if f == 0
            U[i] += dt * Fx
            V[i] += dt * Fy
        else
            # 半隐式 Coriolis：对 (U,V) 做旋转，Stokes–Coriolis 与体力显式
            θ = f * dt
            cθ, sθ = cos(θ), sin(θ)
            Ru = U[i] + dt * (f * Vs + Fx)
            Rv = V[i] + dt * (-f * Us + Fy)
            U[i] = cθ * Ru + sθ * Rv
            V[i] = -sθ * Ru + cθ * Rv
        end
    end
    return U, V
end

"""
由界面有效粘性与表面/底部应力诊断层中心速度剪切。
稳态关系：`νe ∂U/∂z = τ(z)`，其中 `τ(z) = τ_top + Fx*(z - 0) ...` 更一般地由
积分动量得到。此处用中心差分诊断瞬时剪切。
"""
function velocity_shear!(Uz::AbstractVector, Vz::AbstractVector,
                         U::AbstractVector, V::AbstractVector, grid::UniformColumnGrid)
    Nz = grid.Nz
    dz = grid.dz
    @inbounds for i in 1:Nz
        if i == 1
            Uz[i] = (U[2] - U[1]) / dz
            Vz[i] = (V[2] - V[1]) / dz
        elseif i == Nz
            Uz[i] = (U[Nz] - U[Nz-1]) / dz
            Vz[i] = (V[Nz] - V[Nz-1]) / dz
        else
            Uz[i] = (U[i+1] - U[i-1]) / (2dz)
            Vz[i] = (V[i+1] - V[i-1]) / (2dz)
        end
    end
    return Uz, Vz
end
