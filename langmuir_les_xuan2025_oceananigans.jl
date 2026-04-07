# Oceananigans approximation to the LES configuration in
# Xuan & Shen (2025), JFM 1023 A4:
# "Resolvent model-based analyses of coherent structures in Langmuir turbulence"
#
# This script aligns the physical set-up with the paper:
#   * unstratified Langmuir turbulence
#   * domain: 8πH × 4πH × H
#   * deep-water Stokes drift with k0 H = 3.5
#   * turbulent Langmuir number La_t = 0.2 or 0.3
#   * friction Reynolds number Re_tau = 1000
#   * surface wind stress + bottom stress-free
#   * adverse pressure gradient balancing the mean momentum
#   * Dynamic Smagorinsky SGS closure
#
# The paper uses a hybrid pseudo-spectral / finite-difference LES code.
# Oceananigans cannot reproduce that discretization exactly, but the script below
# matches the governing physical configuration as closely as possible.

using Oceananigans
using Oceananigans.Units: minute, minutes, hours
using Oceananigans.Units: GiB
using Printf

# -------------------------
# Case and physical control
# -------------------------

architecture = CPU()  # switch to GPU() on a machine with CUDA support

const H = 25.0                 # m, boundary-layer depth
const La_t = 0.2               # paper cases: 0.2 or 0.3
const Re_tau = 1000.0

# A convenient dimensionalization for the idealized paper set-up.
# Only the non-dimensional groups La_t and Re_tau matter here.
const u_star = 0.01            # m s^-1
const nu = u_star * H / Re_tau # m^2 s^-1
const k0 = 3.5 / H             # m^-1, paper uses k0 H = 3.5
const Us0 = u_star / La_t^2    # m s^-1, from La_t = sqrt(u_star / Us0)

const Lx = 8π * H
const Ly = 4π * H

# Paper resolution: 768 × 256 × 1024 with 256 points in the vertical direction.
# Here x, y are horizontal and z is vertical, so the Oceananigans mapping is
# (Nx, Ny, Nz) = (768, 1024, 256).
#
# This grid is very expensive. Use the exact paper grid only for production runs.
const exact_paper_grid = false
const grid_size = exact_paper_grid ? (768, 1024, 256) : (192, 256, 128)
const Nx, Ny, Nz = grid_size

# Approximate the paper's surface/bottom clustering with a symmetric tanh mesh.
const stretching = 3.5
z_faces(k) = begin
    η = (k - 1) / Nz
    ξ = 2 * η - 1
    0.5 * H * (tanh(stretching * ξ) / tanh(stretching) - 1)
end

grid = RectilinearGrid(architecture;
                       topology = (Periodic, Periodic, Bounded),
                       size = (Nx, Ny, Nz),
                       x = (0, Lx),
                       y = (0, Ly),
                       z = z_faces)

# -------------------------
# Deep-water Stokes drift
# -------------------------

stokes_drift(z) = Us0 * exp(2 * k0 * z)
d_stokes_dz(z, t) = 2 * k0 * Us0 * exp(2 * k0 * z)

# -------------------------
# Boundary conditions
# -------------------------

# Oceananigans uses positive-upward flux conventions. A negative surface flux drives
# a positive x-velocity. The body force sign is chosen to oppose the wind-driven flow.
const tau_w = -u_star^2
const F_u = -u_star^2 / H

u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(tau_w),
                                bottom = FluxBoundaryCondition(0.0))

v_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0.0),
                                bottom = FluxBoundaryCondition(0.0))

forcing_u(x, y, z, t) = F_u

# -------------------------
# Model
# -------------------------

closure = (ScalarDiffusivity(ν = nu), DynamicSmagorinsky())

model = NonhydrostaticModel(; grid,
                            coriolis = nothing,
                            advection = WENO(),
                            timestepper = :RungeKutta3,
                            closure,
                            stokes_drift = UniformStokesDrift(∂z_uˢ = d_stokes_dz),
                            boundary_conditions = (u = u_bcs, v = v_bcs),
                            forcing = (u = forcing_u,))

# Oceananigans prognoses the Lagrangian-mean velocity when Stokes drift is enabled.
# Thus model.velocities.u is already U_L in the notation of the paper.

# -------------------------
# Initial condition
# -------------------------

Ξ(z) = randn() * exp(z / (0.15 * H))
u_i(x, y, z) = 1e-3 * u_star * Ξ(z)
v_i(x, y, z) = 1e-3 * u_star * Ξ(z)
w_i(x, y, z) = 1e-3 * u_star * Ξ(z)

set!(model, u = u_i, v = v_i, w = w_i)

# -------------------------
# Simulation set-up
# -------------------------

simulation = Simulation(model, Δt = 0.25, stop_time = 72hours)

wizard = TimeStepWizard(cfl = 0.8, max_change = 1.1, max_Δt = 0.5minute)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

function progress(simulation)
    u, v, w = simulation.model.velocities

    msg = @sprintf("i: %06d, t: %s, Δt: %s, umax = (%.3e, %.3e, %.3e) m s^-1, wall time: %s",
                   iteration(simulation),
                   prettytime(time(simulation)),
                   prettytime(simulation.Δt),
                   maximum(abs, u),
                   maximum(abs, v),
                   maximum(abs, w),
                   prettytime(simulation.run_wall_time))

    @info msg
    return nothing
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(50))

# -------------------------
# Diagnostics for comparison with the paper
# -------------------------

u_L = model.velocities.u
v = model.velocities.v
w = model.velocities.w

stokes_velocity = Field{Face, Center, Center}(grid)
set!(stokes_velocity, (x, y, z) -> stokes_drift(z))

u_E = Field(u_L - stokes_velocity)

U_L = Field(Average(u_L, dims = (1, 2)))
U_E = Field(Average(u_E, dims = (1, 2)))
uw = Field(Average(@at (Face, Center, Face) u_L * w, dims = (1, 2)))
vw = Field(Average(@at (Center, Face, Face) v * w, dims = (1, 2)))

function eddy_viscosity_field(model)
    if hasproperty(model, :closure_fields) && hasproperty(model.closure_fields, :νₑ)
        return model.closure_fields.νₑ
    else
        return model.diffusivity_fields.νₑ
    end
end

nu_e = eddy_viscosity_field(model)
nu_e_bar = Field(Average(nu_e, dims = (1, 2)))

output_interval = 5minutes

simulation.output_writers[:fields] =
    JLD2OutputWriter(model, (; u_L, u_E, v, w, nu_e),
                     schedule = TimeInterval(output_interval),
                     filename = "langmuir_les_fields.jld2",
                     max_filesize = 8GiB,
                     overwrite_existing = true)

simulation.output_writers[:averages] =
    JLD2OutputWriter(model, (; U_L, U_E, uw, vw, nu_e_bar),
                     schedule = AveragedTimeInterval(output_interval, window = 2minutes),
                     filename = "langmuir_les_averages.jld2",
                     max_filesize = 8GiB,
                     overwrite_existing = true)

simulation.output_writers[:checkpointer] =
    Checkpointer(model,
                 schedule = TimeInterval(12hours),
                 prefix = "langmuir_les_checkpoint",
                 overwrite_existing = true)

run!(simulation)
