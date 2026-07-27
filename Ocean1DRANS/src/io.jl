"""
将稳态解导出为字典（含 Lagrangian 平均流 UL = U + Us）。
"""
function profile_dict(sol::SteadySolution)
    cfg = sol.config
    st = sol.state
    g = cfg.grid
    UL = st.U .+ cfg.stokes.us_c
    VL = st.V .+ cfg.stokes.vs_c
    return Dict(
        "z" => copy(g.zc),
        "zf" => copy(g.zf),
        "U" => copy(st.U),
        "V" => copy(st.V),
        "UL" => UL,
        "VL" => VL,
        "k" => copy(st.k),
        "ell" => copy(st.ℓ),
        "nu_t" => copy(st.νt_c),
        "nu_cl" => copy(st.νcl_c),
        "nu_t_faces" => copy(st.νt_f),
        "Us" => copy(cfg.stokes.us_c),
        "Vs" => copy(cfg.stokes.vs_c),
        "dUsdz" => copy(cfg.stokes.dusdz_c),
        "La_t" => cfg.La_t,
        "k0H" => cfg.k0H,
        "u_star" => cfg.forcing.u★,
        "H" => g.H,
        "nu" => cfg.forcing.ν,
        "alpha_s" => closure_αs(cfg.closure),
        "converged" => sol.converged,
        "residual" => sol.residual,
        "iterations" => sol.iterations,
    )
end

"""
    write_profiles_csv(path, sol)

写出层中心廓线：`z, U, V, UL, VL, Us, Vs, k, nu_t, ell`。

论文 Fig.2(a) 对应的是 **UL**（Lagrangian），不是欧拉 U。
"""
function write_profiles_csv(path::AbstractString, sol::SteadySolution)
    cfg = sol.config
    st = sol.state
    g = cfg.grid
    open(path, "w") do io
        println(io, "z,U,V,UL,VL,Us,Vs,k,nu_t,nu_cl,ell")
        for i in eachindex(g.zc)
            UL = st.U[i] + cfg.stokes.us_c[i]
            VL = st.V[i] + cfg.stokes.vs_c[i]
            @printf(io, "%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e\n",
                    g.zc[i], st.U[i], st.V[i], UL, VL,
                    cfg.stokes.us_c[i], cfg.stokes.vs_c[i],
                    st.k[i], st.νt_c[i], st.νcl_c[i], st.ℓ[i])
        end
    end
    return path
end
