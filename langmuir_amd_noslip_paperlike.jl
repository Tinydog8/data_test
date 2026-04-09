# -*- coding: utf-8 -*-

ENV["CUDA_VISIBLE_DEVICES"] = "2"

using Oceananigans
using Oceananigans.Units: minute, minutes, hours
using Oceananigans.Units: GiB
using DataFrames
using Printf
using CUDA

# ---------------------------------------------------
# Paper-like Langmuir set-up with bottom no-slip wall
# ---------------------------------------------------
#
# This script follows Xuan & Shen (2025) as closely as practical in Oceananigans,
# but changes the bottom boundary from stress-free to a no-slip solid wall.
#
# Main differences relative to the paper:
#   1. SGS closure is AMD instead of dynamic Smagorinsky.
#   2. Bottom boundary is no-slip rather than stress-free.
#   3. The original LES solver/discretization is not reproduced exactly.
#
# Main paper-like choices preserved:
#   * deep-water Stokes drift with k0 * H = 3.5
#   * domain size 8πH × 4πH × H
#   * La_t = 0.2 or 0.3 control
#   * Re_tau = u★ H / ν = 1000
#   * top wind stress, horizontally periodic domain

# -------------------------
# Domain and control params
# -------------------------

const architecture = GPU()

const wavelength = 60.0                    # m
const wavenumber = 2π / wavelength         # m^-1
const H = 3.5 / wavenumber                 # m, from paper k0 * H = 3.5
const Lx = 8π * H
const Ly = 4π * H

const La_t = 0.2                           # paper uses 0.2 and 0.3
const Reτ = 1000.0

# Physical molecular viscosity of seawater. Using this with Re_tau = 1000 implies
# a very small u★. This is mathematically consistent, but the flow may require a
# very long physical spin-up time to reach a turbulent Langmuir state.
const νₘ = 1.05e-6

# Paper-like wave setting via deep-water monochromatic Stokes drift.
const Uˢ₀ = 0.15                           # m s^-1, representative surface Stokes drift
const u★ = La_t^2 * Uˢ₀                    # from La_t = sqrt(u★ / Uˢ₀)
const ν_target = u★ * H / Reτ

# Keep the paper-like non-dimensional control exact while reporting the implied
# Reynolds number under physical seawater viscosity.
const Reτ_physical_ν = u★ * H / νₘ
const Qᵘ = -u★^2

# Use a practical grid that respects the paper's aspect ratio and bottom/surface
# refinement idea without reproducing the exact original LES discretization.
const Nx = 600
const Ny = 300
const Nz = 128
const stretching = 2.5

z_faces(k) = begin
    η = (k - 1) / Nz
    ξ = 2η - 1
    0.5 * H * (tanh(stretching * ξ) / tanh(stretching) - 1)
end

grid = RectilinearGrid(architecture;
                       topology = (Periodic, Periodic, Bounded),
                       size = (Nx, Ny, Nz),
                       x = (0, Lx),
                       y = (0, Ly),
                       z = z_faces)

# -------------------------
# Stokes drift
# -------------------------

uˢ(z) = Uˢ₀ * exp(2 * wavenumber * z)
∂z_uˢ(z, t) = 2 * wavenumber * Uˢ₀ * exp(2 * wavenumber * z)

# -------------------------
# Boundary conditions
# -------------------------

# We interpret "bottom no-slip" as Eulerian no-slip. Because Oceananigans evolves
# the Lagrangian-mean velocity when Stokes drift is present, a strict Eulerian
# no-slip condition corresponds to setting u_L(bottom) = uˢ(bottom).
const u_bottom_value = uˢ(-H)

u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵘ),
                                bottom = ValueBoundaryCondition(u_bottom_value))

v_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0.0),
                                bottom = ValueBoundaryCondition(0.0))

# -------------------------
# Diagnostics metadata
# -------------------------

data = Dict(
    "H" => H,
    "Lx" => Lx,
    "Ly" => Ly,
    "wavelength" => wavelength,
    "k0H" => wavenumber * H,
    "La_t" => La_t,
    "u_star" => u★,
    "surface_stress_Q_u" => Qᵘ,
    "surface_stokes_drift" => Uˢ₀,
    "Re_tau_target" => Reτ,
    "nu_for_target_Re_tau" => ν_target,
    "physical_nu" => νₘ,
    "Re_tau_with_physical_nu" => Reτ_physical_ν,
    "bottom_uL_value" => u_bottom_value,
)

println(DataFrame(data))

# -------------------------
# Model
# -------------------------

coriolis = nothing

# Keep AMD as requested, but add explicit molecular viscosity so that the target
# Re_tau is actually realized in the momentum equation.
closure = (ScalarDiffusivity(ν = ν_target), AnisotropicMinimumDissipation())

model = NonhydrostaticModel(; grid, coriolis,
                            advection = WENO(),
                            timestepper = :RungeKutta3,
                            closure,
                            stokes_drift = UniformStokesDrift(∂z_uˢ = ∂z_uˢ),
                            boundary_conditions = (u = u_bcs, v = v_bcs))

# -------------------------
# Initial conditions
# -------------------------

Ξ(z) = randn() * exp(z / (0.12 * H))

uᵢ(x, y, z) = u_bottom_value + 1e-2 * u★ * Ξ(z)
vᵢ(x, y, z) = 1e-2 * u★ * Ξ(z)
wᵢ(x, y, z) = 1e-3 * u★ * Ξ(z)

set!(model, u = uᵢ, v = vᵢ, w = wᵢ)

# -------------------------
# Simulation control
# -------------------------

simulation = Simulation(model, Δt = 0.25, stop_time = 72hours)

wizard = TimeStepWizard(cfl = 0.5, max_change = 1.05, max_Δt = 20.0)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(1))

function progress(simulation)
    u, v, w = simulation.model.velocities

    msg = @sprintf("i: %06d, t: %s, Δt: %s, umax = (%.5e, %.5e, %.5e) ms⁻¹, wall time: %s",
                   iteration(simulation),
                   prettytime(time(simulation)),
                   prettytime(simulation.Δt),
                   maximum(abs, u), maximum(abs, v), maximum(abs, w),
                   prettytime(simulation.run_wall_time))

    @info msg
    return nothing
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(50))

# -------------------------
# Diagnostics
# -------------------------

u_L = model.velocities.u
v = model.velocities.v
w = model.velocities.w

stokes_velocity = Field{Face, Center, Center}(grid)
set!(stokes_velocity, (x, y, z) -> uˢ(z))

u_E = Field(u_L - stokes_velocity)

U_L = Field(Average(u_L, dims = (1, 2)))
U_E = Field(Average(u_E, dims = (1, 2)))
uw = Field(Average(@at (Face, Center, Face) u_L * w, dims = (1, 2)))
vw = Field(Average(@at (Center, Face, Face) v * w, dims = (1, 2)))

function eddy_viscosity_field(model)
    if hasproperty(model, :closure_fields) && hasproperty(model.closure_fields, :νₑ)
        return model.closure_fields.νₑ
    elseif hasproperty(model, :diffusivity_fields) && hasproperty(model.diffusivity_fields, :νₑ)
        return model.diffusivity_fields.νₑ
    else
        return last(model.closure_fields).νₑ
    end
end

nu_e = eddy_viscosity_field(model)
nu_e_bar = Field(Average(nu_e, dims = (1, 2)))

dudz = Field(Average(@at (Center, Center, Face) ∂z(u_L), dims = (1, 2)))

output_interval = 5minutes

simulation.output_writers[:averages] =
    JLD2Writer(model, (; U_L, U_E, uw, vw, nu_e_bar, dudz),
               schedule = AveragedTimeInterval(output_interval, window = 2minutes),
               filename = "langmuir_noslip_amd_averages.jld2",
               overwrite_existing = true)

simulation.output_writers[:vertical_velocity] =
    JLD2Writer(model, (; w),
               schedule = TimeInterval(output_interval),
               filename = "langmuir_noslip_amd_w.jld2",
               overwrite_existing = true)

simulation.output_writers[:checkpointer] =
    Checkpointer(model,
                 schedule = TimeInterval(12hours),
                 prefix = "langmuir_noslip_amd_checkpoint",
                 overwrite_existing = true)

run!(simulation)
