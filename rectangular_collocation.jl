using LinearAlgebra

const RECT_COLLOC_STAMP = "rectcolloc-2026-04-03"

struct RectGrid
    N::Int
    H::Float64
    y_int::Vector{Float64}
    y_v::Vector{Float64}
    y_w::Vector{Float64}
    w_y_int::Vector{Float64}
    Dint::Matrix{Float64}
    Dv::Matrix{Float64}
    D2v::Matrix{Float64}
    Dw::Matrix{Float64}
    D2w::Matrix{Float64}
    Pv::Matrix{Float64}
    Pw::Matrix{Float64}
end

function chebyshev_second_kind_nodes(n::Integer)
    n >= 2 || error("second-kind Chebyshev nodes require n >= 2")
    j = 0:(n - 1)
    return -cos.(pi .* j ./ (n - 1))
end

function chebyshev_first_kind_nodes(n::Integer)
    n >= 1 || error("first-kind Chebyshev nodes require n >= 1")
    j = 0:(n - 1)
    return -cos.((j .+ 0.5) .* pi ./ n)
end

function barycentric_weights(x::AbstractVector{<:Real})
    n = length(x)
    λ = ones(Float64, n)
    for j in 1:n
        xj = Float64(x[j])
        prod = 1.0
        for k in 1:n
            k == j && continue
            prod *= (xj - Float64(x[k]))
        end
        λ[j] = inv(prod)
    end
    λ ./= maximum(abs, λ)
    return λ
end

function differentiation_matrix(x::AbstractVector{<:Real}, λ::AbstractVector{<:Real})
    n = length(x)
    length(λ) == n || error("x and λ length mismatch")
    D = zeros(Float64, n, n)
    for i in 1:n
        xi = Float64(x[i])
        λi = Float64(λ[i])
        rowsum = 0.0
        for j in 1:n
            i == j && continue
            val = Float64(λ[j]) / λi / (xi - Float64(x[j]))
            D[i, j] = val
            rowsum += val
        end
        D[i, i] = -rowsum
    end
    return D
end

function resampling_matrix(
    xsrc::AbstractVector{<:Real},
    xdst::AbstractVector{<:Real},
    λsrc::AbstractVector{<:Real};
    atol::Float64 = 1e-14,
)
    nsrc = length(xsrc)
    length(λsrc) == nsrc || error("xsrc and λsrc length mismatch")
    P = zeros(Float64, length(xdst), nsrc)
    for i in eachindex(xdst)
        x = Float64(xdst[i])
        assigned = false
        for j in 1:nsrc
            if abs(x - Float64(xsrc[j])) <= atol
                P[i, j] = 1.0
                assigned = true
                break
            end
        end
        assigned && continue
        denom = 0.0
        for j in 1:nsrc
            denom += Float64(λsrc[j]) / (x - Float64(xsrc[j]))
        end
        for j in 1:nsrc
            P[i, j] = (Float64(λsrc[j]) / (x - Float64(xsrc[j]))) / denom
        end
    end
    return P
end

function gauss_legendre(n::Integer)
    n >= 1 || error("Gauss-Legendre order must be positive")
    if n == 1
        return [0.0], [2.0]
    end
    β = [k / sqrt(4.0 * k^2 - 1.0) for k in 1:(n - 1)]
    J = SymTridiagonal(zeros(Float64, n), β)
    eig = eigen(J)
    x = Vector{Float64}(eig.values)
    w = 2.0 .* abs2.(eig.vectors[1, :])
    return x, w
end

function interpolatory_quadrature_weights(
    xnodes::AbstractVector{<:Real},
    λnodes::AbstractVector{<:Real},
    nquad::Integer,
)
    xq, wq = gauss_legendre(nquad)
    Pq = resampling_matrix(xnodes, xq, λnodes)
    return vec(Pq' * wq)
end

scale_from_reference(x::AbstractVector{<:Real}, H::Real) = 0.5 .* Float64(H) .* (Float64.(x) .- 1.0)

function build_rect_grid(N::Integer, H::Real)
    N >= 8 || error("N should be at least 8 for the fourth-/second-order system")
    Hf = Float64(H)

    Nv = N + 4
    Nw = N + 2

    x_v = chebyshev_second_kind_nodes(Nv)
    x_w = chebyshev_second_kind_nodes(Nw)
    x_int = chebyshev_first_kind_nodes(N)

    λ_v = barycentric_weights(x_v)
    λ_w = barycentric_weights(x_w)
    λ_int = barycentric_weights(x_int)

    Dv_ref = differentiation_matrix(x_v, λ_v)
    Dw_ref = differentiation_matrix(x_w, λ_w)
    Dint_ref = differentiation_matrix(x_int, λ_int)

    scale = 2.0 / Hf
    Dv = scale .* Dv_ref
    Dw = scale .* Dw_ref
    Dint = scale .* Dint_ref

    D2v = Dv * Dv
    D2w = Dw * Dw

    Pv = resampling_matrix(x_v, x_int, λ_v)
    Pw = resampling_matrix(x_w, x_int, λ_w)

    # Use an interpolatory quadrature on the interior grid to approximate (2.21).
    w_ref = interpolatory_quadrature_weights(x_int, λ_int, max(2N, 256))
    w_y = 0.5 * Hf .* w_ref

    return RectGrid(
        N,
        Hf,
        scale_from_reference(x_int, Hf),
        scale_from_reference(x_v, Hf),
        scale_from_reference(x_w, Hf),
        w_y,
        Dint,
        Dv,
        D2v,
        Dw,
        D2w,
        Pv,
        Pw,
    )
end

function transfer_gain_rect(
    kx::Real,
    kz::Real,
    omega::Real,
    g::RectGrid,
    Uv::AbstractVector,
    Usv::AbstractVector,
    nuTv::AbstractVector,
    dUv::AbstractVector,
    d2Uv::AbstractVector,
    dUsv::AbstractVector,
    dnuTv::AbstractVector,
    d2nuTv::AbstractVector,
    Uw::AbstractVector,
    Usw::AbstractVector,
    nuTw::AbstractVector,
    dUw::AbstractVector,
    d2Uw::AbstractVector,
    dUsw::AbstractVector,
    dnuTw::AbstractVector,
    d2nuTw::AbstractVector,
)
    N = g.N
    Nv = length(g.y_v)
    Nw = length(g.y_w)

    for (name, arr, expected) in (
        ("Uv", Uv, Nv),
        ("Usv", Usv, Nv),
        ("nuTv", nuTv, Nv),
        ("dUv", dUv, Nv),
        ("d2Uv", d2Uv, Nv),
        ("dUsv", dUsv, Nv),
        ("dnuTv", dnuTv, Nv),
        ("d2nuTv", d2nuTv, Nv),
        ("Uw", Uw, Nw),
        ("Usw", Usw, Nw),
        ("nuTw", nuTw, Nw),
        ("dUw", dUw, Nw),
        ("d2Uw", d2Uw, Nw),
        ("dUsw", dUsw, Nw),
        ("dnuTw", dnuTw, Nw),
        ("d2nuTw", d2nuTw, Nw),
    )
        length(arr) == expected || error("$name has length $(length(arr)), expected $expected")
    end

    k2 = Float64(kx)^2 + Float64(kz)^2
    k2 > 0 || error("k_x = k_z = 0 makes C singular; choose at least one non-zero horizontal wavenumber")
    ikx = im * Float64(kx)
    ikz = im * Float64(kz)

    Iv = Matrix{Float64}(I, Nv, Nv)
    Iw = Matrix{Float64}(I, Nw, Nw)
    Ii = Matrix{Float64}(I, N, N)

    Δv = g.D2v - k2 .* Iv
    Δw = g.D2w - k2 .* Iw

    ULv = Float64.(Uv) .+ Float64.(Usv)
    ULw = Float64.(Uw) .+ Float64.(Usw)

    # Discrete version of (2.13) and (2.14), with the Appendix-B rectangular resampling.
    LOS = g.Pv * (
        -ikx .* Diagonal(ULv) * Δv +
        ikx .* Diagonal(Float64.(d2Uv)) +
        Diagonal(Float64.(nuTv)) * (Δv * Δv) +
        2.0 .* Diagonal(Float64.(dnuTv)) * g.Dv * Δv +
        2.0 .* Diagonal(Float64.(d2nuTv)) * (g.D2v + k2 .* Iv)
    )
    LSq = g.Pw * (
        -ikx .* Diagonal(ULw) +
        Diagonal(Float64.(nuTw)) * Δw +
        Diagonal(Float64.(dnuTw)) * g.Dw
    )

    Usprime_int = g.Pv * Float64.(dUsv)
    Uprime_int = g.Pw * Float64.(dUw)

    L12 = -ikz .* Diagonal(Usprime_int) * g.Pw
    L21 = -ikz .* Diagonal(Uprime_int) * g.Pv

    E11 = g.Pv * Δv
    E22 = g.Pw

    nunk = Nv + Nw
    A = zeros(ComplexF64, nunk, nunk)
    A[1:N, 1:Nv] = im * Float64(omega) .* E11 - LOS
    A[1:N, (Nv + 1):end] = -L12
    A[(N + 1):(2N), 1:Nv] = -L21
    A[(N + 1):(2N), (Nv + 1):end] = im * Float64(omega) .* E22 - LSq

    # Boundary conditions: v = D^2 v = 0, Dω_y = 0 at y = 0 and y = -H.
    A[(2N + 1), 1:Nv] = Iv[1, :]
    A[(2N + 2), 1:Nv] = Iv[end, :]
    A[(2N + 3), 1:Nv] = g.D2v[1, :]
    A[(2N + 4), 1:Nv] = g.D2v[end, :]
    A[(2N + 5), (Nv + 1):end] = g.Dw[1, :]
    A[(2N + 6), (Nv + 1):end] = g.Dw[end, :]

    Dv_int = g.Pv * g.Dv

    C = zeros(ComplexF64, 3N, nunk)
    C[1:N, 1:Nv] = (ikx / k2) .* Dv_int
    C[1:N, (Nv + 1):end] = (-ikz / k2) .* g.Pw
    C[(N + 1):(2N), 1:Nv] = g.Pv
    C[(2N + 1):(3N), 1:Nv] = (ikz / k2) .* Dv_int
    C[(2N + 1):(3N), (Nv + 1):end] = (ikx / k2) .* g.Pw

    # Force components are discretised on the interior first-kind grid.
    B = zeros(ComplexF64, nunk, 3N)
    B[1:N, 1:N] = (-ikx) .* g.Dint
    B[1:N, (N + 1):(2N)] = (-k2) .* Ii
    B[1:N, (2N + 1):(3N)] = (-ikz) .* g.Dint
    B[(N + 1):(2N), 1:N] = ikz .* Ii
    B[(N + 1):(2N), (2N + 1):(3N)] = (-ikx) .* Ii

    T = C * (A \ B)
    svals = svdvals(T)
    return (
        T = T,
        A = A,
        B = B,
        C = C,
        E11 = E11,
        E22 = E22,
        G = abs2(svals[1]),
        sigma = svals[1],
    )
end
