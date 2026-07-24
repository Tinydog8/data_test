"""
将稳态解导出为字典（便于绘图 / 与 resolvent 代码对接）。
"""
function profile_dict(sol::SteadySolution)
    cfg = sol.config
    st = sol.state
    g = cfg.grid
    return Dict(
        "z" => copy(g.zc),
        "zf" => copy(g.zf),
        "U" => copy(st.U),
        "V" => copy(st.V),
        "k" => copy(st.k),
        "ell" => copy(st.ℓ),
        "nu_t" => copy(st.νt_c),
        "nu_t_faces" => copy(st.νt_f),
        "Us" => copy(cfg.stokes.us_c),
        "Vs" => copy(cfg.stokes.vs_c),
        "dUsdz" => copy(cfg.stokes.dusdz_c),
        "La_t" => cfg.La_t,
        "k0H" => cfg.k0H,
        "u_star" => cfg.forcing.u★,
        "H" => g.H,
        "nu" => cfg.forcing.ν,
        "converged" => sol.converged,
        "residual" => sol.residual,
        "iterations" => sol.iterations,
    )
end

"""
    write_profiles_csv(path, sol)

写出层中心廓线 CSV：`z, U, V, Us, Vs, k, nu_t, ell`。
"""
function write_profiles_csv(path::AbstractString, sol::SteadySolution)
    cfg = sol.config
    st = sol.state
    g = cfg.grid
    open(path, "w") do io
        println(io, "z,U,V,Us,Vs,k,nu_t,ell")
        for i in eachindex(g.zc)
            @printf(io, "%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e\n",
                    g.zc[i], st.U[i], st.V[i],
                    cfg.stokes.us_c[i], cfg.stokes.vs_c[i],
                    st.k[i], st.νt_c[i], st.ℓ[i])
        end
    end
    return path
end
