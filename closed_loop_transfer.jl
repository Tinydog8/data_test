using LinearAlgebra
using Logging

# =========================
# Chebyshev differentiation
# =========================
function chebydif(N::Int)
    x = cos.(pi .* collect(0:N) ./ N)              # x goes from 1 -> -1
    D = zeros(Float64, N + 1, N + 1)
    for i in 0:N
        ci = (i == 0 || i == N) ? 2.0 : 1.0
        xi = x[i + 1]
        for j in 0:N
            cj = (j == 0 || j == N) ? 2.0 : 1.0
            if i == j && i == 0
                D[i + 1, j + 1] = (2.0 * N^2 + 1.0) / 6.0
            elseif i == j && i == N
                D[i + 1, j + 1] = -(2.0 * N^2 + 1.0) / 6.0
            elseif i == j
                D[i + 1, j + 1] = -xi / (2.0 * (1.0 - xi^2))
            else
                D[i + 1, j + 1] = (ci / cj) * (-1)^(i + j) / (xi - x[j + 1])
            end
        end
    end
    return x, D
end

# =========================
# Base flows
# =========================
U_couette(z) = z
U_poiseuille(z) = 1.0 .- z.^2
U_prime_couette(z) = ones(length(z))
U_prime_poiseuille(z) = -2.0 .* z

# =========================
# Stokes drift profiles
# =========================
function stokes_drift_profile(z, u_s0_star, k_w_star; depth_scale=1.0)
    # z ∈ [-1, 1], z=1 is surface, z=-1 is bottom
    # Use decay with distance from surface (1 - z)
    return u_s0_star .* exp.(-2.0 .* k_w_star .* (1.0 .- z) ./ depth_scale)
end

# =========================
# Matrix builders (IMPORTANT: keep z and D consistent)
# =========================
function build_matrices(N, k; flow_type::String="poiseuille")
    z_raw, D_raw = chebydif(N)                 # z_raw: 1 -> -1
    n = N + 1

    # Reorder to z increasing: -1 -> 1, and reorder D consistently
    p = n:-1:1
    z = z_raw[p]
    D1 = D_raw[p, p]

    I = Matrix{Float64}(I, n, n)
    D2 = D1 * D1
    Δ = D2 .- k^2 .* I
    Δ2 = Δ * Δ

    U′ = flow_type == "couette" ? U_prime_couette(z) : U_prime_poiseuille(z)
    return z, D1, Δ, Δ2, U′
end

# ==========================================
# Closed-loop transfer: block system solution
# ==========================================
"""
    G_closed_frequency(Δ, Δ2, D1, U′, u_s, La, k, Re, ω) -> G

Compute closed-loop transfer matrix `G` mapping stacked forcing `[f_x; f_y; f_z]`
to wall-normal velocity `w` at frequency ω (s = iω), i.e.

    ŵ = G( iω ) * [ f̂_x; f̂_y; f̂_z ] .

Discretization uses Chebyshev collocation with tau row replacement for BC:
- u(±1)=0
- w(±1)=0, Dw(±1)=0
"""
function G_closed_frequency(Δ, Δ2, D1, U′, u_s, La, k, Re, ω)
    s = im * ω
    n = size(Δ, 1)

    I = Matrix{ComplexF64}(I, n, n)
    Δc = ComplexF64.(Δ)
    Δ2c = ComplexF64.(Δ2)
    D1c = ComplexF64.(D1)

    # u_s as vector on grid
    u_s_vec = u_s isa Number ? fill(Float64(u_s), n) : collect(Float64.(u_s))
    u_s′ = D1 * u_s_vec                               # u_s'(z) on grid
    Us′diag = Diagonal(ComplexF64.(u_s′))
    U′diag = Diagonal(ComplexF64.(U′))

    # Block operator for [u; w]
    Auu = s .* I .- Δc ./ Re
    Auw = U′diag
    Awu = -(1.0 / La) * k^2 .* Us′diag
    Aww = s .* Δc .- Δ2c ./ Re

    A = [Auu  Auw;
         Awu  Aww]

    # RHS mapping for stacked forcing [fx; fy; fz]
    Z = zeros(ComplexF64, n, n)
    Fx = I
    Fy = (-im * k) .* D1c
    Fz = (-k^2) .* I
    RHS = [Fx  Z   Z;
           Z   Fy  Fz]    # 2n x 3n

    # Enforce BC via row replacement (tau)
    # u(±1)=0 -> rows 1 and n of the u-block equation
    A[1, :] .= 0
    A[1, 1:n] .= I[1, :]
    RHS[1, :] .= 0

    A[n, :] .= 0
    A[n, 1:n] .= I[n, :]
    RHS[n, :] .= 0

    # w(±1)=0, Dw(±1)=0 -> rows in w-block (offset by n)
    # w(-1)=0 at row n+1, w(1)=0 at row 2n
    A[n + 1, :] .= 0
    A[n + 1, n + 1:2n] .= I[1, :]
    RHS[n + 1, :] .= 0

    A[2n, :] .= 0
    A[2n, n + 1:2n] .= I[n, :]
    RHS[2n, :] .= 0

    # Dw(-1)=0 at row n+2 uses D1 first-row; Dw(1)=0 at row 2n-1 uses D1 last-row
    A[n + 2, :] .= 0
    A[n + 2, n + 1:2n] .= D1c[1, :]
    RHS[n + 2, :] .= 0

    A[2n - 1, :] .= 0
    A[2n - 1, n + 1:2n] .= D1c[n, :]
    RHS[2n - 1, :] .= 0

    # Solve once for the whole transfer matrix
    F = lu(A) \ RHS
    Wmap = F[n + 1:2n, :]     # extract w rows
    return Wmap               # n x (3n)
end

# =========================
# H∞ norm (frequency sweep)
# =========================
function hinf_norm(G_func, ω_range, args...)
    maxσ = 0.0
    bestω = first(ω_range)
    failed = 0

    for ω in ω_range
        try
            G = G_func(args..., ω)
            if any(!isfinite, G)
                failed += 1
                continue
            end
            σ = maximum(svdvals(Matrix(G)))
            if isfinite(σ) && σ > maxσ
                maxσ = σ
                bestω = ω
            end
        catch
            failed += 1
            continue
        end
    end

    if failed > 0
        @warn "Failed at $failed / $(length(ω_range)) frequencies."
    end
    return maxσ, bestω
end

# ==================================
# Convenience: sweep k for closed-loop
# ==================================
function compute_hinf_vs_k(k_values; flow_type="couette", N=200, Re=1.0, La=0.1,
                           u_s_params=(1.0, 2pi / 2.4), ω_range=range(-5.0, 5.0, length=201))
    norms = Float64[]
    bestωs = Float64[]

    for k in k_values
        z, D1, Δ, Δ2, U′ = build_matrices(N, k; flow_type=flow_type)

        u_s = if u_s_params isa Number
            u_s_params
        elseif u_s_params isa Tuple
            u_s0_star, k_w_star = u_s_params
            stokes_drift_profile(z, u_s0_star, k_w_star)
        elseif u_s_params isa Function
            u_s_params(z)
        else
            error("u_s_params must be Number, Tuple(u_s0_star,k_w_star), or Function")
        end

        σ, ω⋆ = hinf_norm((Δ, Δ2, D1, U′, u_s, La, k, Re, ω) -> G_closed_frequency(Δ, Δ2, D1, U′, u_s, La, k, Re, ω),
                         ω_range, Δ, Δ2, D1, U′, u_s, La, k, Re)

        push!(norms, σ)
        push!(bestωs, ω⋆)
    end
    return k_values, norms, bestωs
end

