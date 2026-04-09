# -*- coding: utf-8 -*-

ENV["CUDA_VISIBLE_DEVICES"] = "2"

using Oceananigans
using Oceananigans.Units: minute, minutes, hours
using Oceananigans.Units: GiB
using Oceananigans.BuoyancyModels: g_Earth
using DataFrames
using Printf

# ---------------------------------
# Domain and low-Re control
# ---------------------------------

const H = 25.0 # m, water depth
const architecture = GPU()

# Lower the Reynolds number by keeping the molecular viscosity at a physical
# seawater value and adjusting the friction velocity u★ so that
# Reτ = u★ H / νₘ matches the target value.
const Reτ_target = 250.0

grid = RectilinearGrid(architecture, size = (600, 600, 100), extent = (600, 600, H))

# ---------------------------------
# Stokes drift
# ---------------------------------

const amplitude = 1.0
const wavelength = 60.0
const wavenumber = 2π / wavelength
const frequency = sqrt(g_Earth * wavenumber)

const vertical_scale = wavelength / 4π
const Uˢ = amplitude^2 * wavenumber * frequency

# The original script intended to use a finite-depth Stokes profile. The
# denominator needs parentheses to avoid Julia's left-to-right precedence rules.
uˢ(z) = Uˢ * cosh(2 * wavenumber * (z + H)) / (2 * (sinh(wavenumber * H))^2)
∂z_uˢ(z, t) = Uˢ * wavenumber * sinh(2 * wavenumber * (z + H)) / (sinh(wavenumber * H))^2

# ---------------------------------
# Wind stress, pressure gradient, and bottom drag
# ---------------------------------

const ρₒ = 1026.0

# Approximate physical kinematic viscosity of seawater near room temperature.
const νₘ = 1.05e-6

# Set the surface stress through the target friction Reynolds number.
const u★ = Reτ_target * νₘ / H
const Qᵘ = -u★^2
const Lat = sqrt(u★ / Uˢ)

const up = 0.01
const PGF = up^2 / H
Fx(x, y, z, t) = PGF

z₀ = 0.01
κ = 0.4
z₁ = -0.5 * znodes(Center, grid)[grid.Nz]
cᴰ = (κ / log(z₁ / z₀))^2

@inline drag_u(x, y, t, u, v, p) = -p.cᴰ * sqrt(u^2 + v^2) * u
@inline drag_v(x, y, t, u, v, p) = -p.cᴰ * sqrt(u^2 + v^2) * v

drag_bc_u = FluxBoundaryCondition(drag_u, field_dependencies = (:u, :v), parameters = (; cᴰ))
drag_bc_v = FluxBoundaryCondition(drag_v, field_dependencies = (:u, :v), parameters = (; cᴰ))

u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵘ), bottom = drag_bc_u)
v_bcs = FieldBoundaryConditions(bottom = drag_bc_v)

const θ = 1 / sqrt(1 + up^2 / u★^2)
const uᵇ = sqrt(u★^2 + up^2)
const u_b_steady = uᵇ / sqrt(cᴰ)
const λ_H = wavelength / H
const kH = wavenumber * H

data = Dict(
    "Lat" => Lat,
    "Θ" => θ,
    "kH" => kH,
    "Wave amplitude" => amplitude,
    "Wave length" => wavelength,
    "PGF velocity" => up,
    "PGF" => PGF,
    "u_star" => u★,
    "Re_tau_target" => Reτ_target,
    "molecular_viscosity" => νₘ,
    "surface_stress_Q_u" => Qᵘ,
)

df_data = DataFrame(data)
println(df_data)

# ---------------------------------
# Tracers and buoyancy
# ---------------------------------

Qʰ = 0.0
cᴾ = 4000.0
Qᵀ = Qʰ / (ρₒ * cᴾ)
dTdz = 0.0

T_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵀ),
                                bottom = FluxBoundaryCondition(0.0))
S_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0.0),
                                bottom = FluxBoundaryCondition(0.0))

# ---------------------------------
# Model
# ---------------------------------

coriolis = nothing

# Keep AMD as requested. Many Oceananigans versions assume seawater molecular
# viscosity / diffusivity by default for AMD, while only some versions expose
# ν and κ as constructor keywords. Using the bare constructor is therefore the
# most version-compatible way to keep physical molecular properties here.
closure = AnisotropicMinimumDissipation()

model = NonhydrostaticModel(; grid, coriolis,
                            advection = WENO(),
                            timestepper = :RungeKutta3,
                            tracers = (:T, :S),
                            buoyancy = SeawaterBuoyancy(),
                            closure,
                            stokes_drift = UniformStokesDrift(∂z_uˢ = ∂z_uˢ),
                            boundary_conditions = (u = u_bcs, v = v_bcs, T = T_bcs, S = S_bcs),
                            forcing = (u = Fx,))

# ---------------------------------
# Initial conditions
# ---------------------------------

Ξ(z) = randn() * exp(z / 4)

uᵢ(x, y, z) = 1e-2 * u★ * Ξ(z)
vᵢ(x, y, z) = 1e-2 * u★ * Ξ(z)
wᵢ(x, y, z) = 1e-4 * u★ * Ξ(z)
Tᵢ(x, y, z) = dTdz * z + 1e-8 * Ξ(z) + 290
Sᵢ = 35.0

set!(model, u = uᵢ, v = vᵢ, w = wᵢ, T = Tᵢ, S = Sᵢ)

# ---------------------------------
# Simulation control
# ---------------------------------

simulation = Simulation(model, Δt = 0.4, stop_time = 24hours)

wizard = TimeStepWizard(cfl = 0.8, max_change = 1.1, max_Δt = 0.1minute)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

function progress(simulation)
    u, v, w = simulation.model.velocities

    msg = @sprintf("i: %04d, t: %s, Δt: %s, umax = (%.5e, %.5e, %.5e) ms⁻¹, wall time: %s",
                   iteration(simulation),
                   prettytime(time(simulation)),
                   prettytime(simulation.Δt),
                   maximum(abs, u), maximum(abs, v), maximum(abs, w),
                   prettytime(simulation.run_wall_time))

    @info msg
    return nothing
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(50))

# ---------------------------------
# Diagnostic fields
# ---------------------------------

Stokes_Field = Field{Center, Center, Face}(grid)
Stokes_velocity_field = Field{Center, Center, Center}(grid)
Stokes_center_Field = Field{Center, Center, Center}(grid)

set!(Stokes_Field, (x, y, z) -> Uˢ * wavenumber * sinh(2 * wavenumber * (z + H)) / (sinh(wavenumber * H))^2)
set!(Stokes_velocity_field, (x, y, z) -> Uˢ * cosh(2 * wavenumber * (z + H)) / (2 * (sinh(wavenumber * H))^2))
set!(Stokes_center_Field, (x, y, z) -> Uˢ * wavenumber * sinh(2 * wavenumber * (z + H)) / (sinh(wavenumber * H))^2)

using Oceananigans.BoundaryConditions: getbc
using Oceananigans: fields

@inline kernel_getbc(i, j, k, grid, boundary_condition, clock, fields) =
    getbc(boundary_condition, i, j, grid, clock, fields)

grid = model.grid
clock = model.clock
model_fields = merge(fields(model), model.auxiliary_fields)
u, v, w = model.velocities
T, S = model.tracers

u_bottom_bc = u.boundary_conditions.bottom
v_bottom_bc = v.boundary_conditions.bottom
u_top_bc = u.boundary_conditions.top
v_top_bc = v.boundary_conditions.top
T_top_bc = T.boundary_conditions.top

u_bottom_bc_op = KernelFunctionOperation{Face, Center, Nothing}(kernel_getbc, grid;
                                                                computed_dependencies = (u_bottom_bc, clock, model_fields))
v_bottom_bc_op = KernelFunctionOperation{Center, Face, Nothing}(kernel_getbc, grid;
                                                                computed_dependencies = (v_bottom_bc, clock, model_fields))
u_top_bc_op = KernelFunctionOperation{Face, Center, Nothing}(kernel_getbc, grid;
                                                             computed_dependencies = (u_top_bc, clock, model_fields))
v_top_bc_op = KernelFunctionOperation{Center, Face, Nothing}(kernel_getbc, grid;
                                                             computed_dependencies = (v_top_bc, clock, model_fields))
T_top_bc_op = KernelFunctionOperation{Center, Center, Nothing}(kernel_getbc, grid;
                                                               computed_dependencies = (T_top_bc, clock, model_fields))

top_ubc_field = Field(u_top_bc_op)
top_vbc_field = Field(v_top_bc_op)
top_Tbc_field = Field(T_top_bc_op)
Bottom_ubc_field = Field(u_bottom_bc_op)
Bottom_vbc_field = Field(v_bottom_bc_op)

output_interval = 0.2minute

u_L = model.velocities.u
u_E = Field(u_L - Stokes_velocity_field)

u = model.velocities.u
v = model.velocities.v
w = model.velocities.w
p = model.pressures.pNHS + model.pressures.pHY′
pNHS = model.pressures.pNHS
pHY = model.pressures.pHY′
T = model.tracers.T
S = model.tracers.S
νₑ = model.diffusivity_fields.νₑ
κₑT = model.diffusivity_fields.κₑ.T

U = Field(Average(u, dims = (1, 2)))
V = Field(Average(v, dims = (1, 2)))
W = Field(Average(w, dims = (1, 2)))
P = Field(Average(p, dims = (1, 2)))
Tₐ = Field(Average(T, dims = (1, 2)))

uw_f_op = @at (Face, Center, Face) (u - U) * w
uw_f = Average(uw_f_op, dims = (1, 2))

vw_f_op = @at (Center, Face, Face) (v - V) * w
vw_f = Average(vw_f_op, dims = (1, 2))

wc_f_op = @at (Center, Center, Face) (T - Tₐ) * w
wc_f = Average(wc_f_op, dims = (1, 2))

sgs_uw_op = @at (Face, Center, Face) νₑ * (∂z(u) + ∂x(w))
sgs_uw = Average(sgs_uw_op, dims = (1, 2))

sgs_vw_op = @at (Center, Face, Face) νₑ * (∂z(v) + ∂y(w))
sgs_vw = Average(sgs_vw_op, dims = (1, 2))

sgs_ww_op = @at (Center, Center, Face) 2νₑ * ∂z(w)
sgs_ww = Average(sgs_ww_op, dims = (1, 2))

sgs_wc_op = @at (Center, Center, Face) κₑT * ∂z(T)
sgs_wc = Average(sgs_wc_op, dims = (1, 2))

u_top = Field(top_ubc_field)
v_top = Field(top_vbc_field)
u_bottom = Field(Bottom_ubc_field)
v_bottom = Field(Bottom_vbc_field)
T_top = Field(top_Tbc_field)

S₁₁_op = @at (Center, Center, Face) ∂x(u)
S₂₂_op = @at (Center, Center, Face) ∂y(v)
S₃₃_op = @at (Center, Center, Face) ∂z(w)
S₁₂_op = @at (Center, Center, Face) ∂y(u) + ∂x(v)
S₁₃_op = @at (Center, Center, Face) ∂z(u) + ∂x(w)
S₂₃_op = @at (Center, Center, Face) ∂y(w) + ∂z(v)

τ₁₁ = @at (Center, Center, Face) Field(νₑ * S₁₁_op)
τ₂₂ = @at (Center, Center, Face) Field(νₑ * S₂₂_op)
τ₃₃ = @at (Center, Center, Face) Field(νₑ * S₃₃_op)
τ₁₂ = @at (Center, Center, Face) Field(νₑ * S₁₂_op)
τ₁₃ = @at (Center, Center, Face) Field(νₑ * S₁₃_op)
τ₂₃ = @at (Center, Center, Face) Field(νₑ * S₂₃_op)

τ₁₁A = Field(Average(τ₁₁, dims = (1, 2)))
τ₂₂A = Field(Average(τ₂₂, dims = (1, 2)))
τ₃₃A = Field(Average(τ₃₃, dims = (1, 2)))
τ₁₂A = Field(Average(τ₁₂, dims = (1, 2)))
τ₁₃A = Field(Average(τ₁₃, dims = (1, 2)))
τ₂₃A = Field(Average(τ₂₃, dims = (1, 2)))

pᶠwᶠ = @at (Center, Center, Face) (p - P) * (w - W)
uᶠuᶠwᶠ = @at (Face, Center, Face) (u - U)^2 * w
vᶠvᶠwᶠ = @at (Center, Face, Face) (v - V)^2 * w
wᶠwᶠwᶠ = @at (Center, Face, Face) w^3

dudz_op = @at (Nothing, Nothing, Face) ∂z(U)
dvdz_op = @at (Nothing, Nothing, Face) ∂z(V)
dwdz_op = @at (Nothing, Nothing, Face) ∂z(W)

shear_budget = @at (Center, Center, Center) uw_f_op * dudz_op

dutdz_op = @at (Center, Center, Face) ∂z(uᶠuᶠwᶠ)
dvtdz_op = @at (Center, Center, Face) ∂z(vᶠvᶠwᶠ)
dwtdz_op = @at (Center, Center, Face) ∂z(wᶠwᶠwᶠ)

τ₁₃ᶠ = @at (Center, Center, Face) (τ₁₃ - τ₁₃A)
τ₂₃ᶠ = @at (Center, Center, Face) (τ₂₃ - τ₂₃A)
τ₃₃ᶠ = @at (Center, Center, Face) (τ₃₃ - τ₃₃A)

uᶠτ₁₃ᶠ = @at (Center, Center, Face) τ₁₃ᶠ * (u - U)
vᶠτ₂₃ᶠ = @at (Center, Center, Face) τ₂₃ᶠ * (v - V)
wᶠτ₃₃ᶠ = @at (Center, Center, Face) τ₃₃ᶠ * w

duτdz_op = @at (Center, Center, Face) ∂z(uᶠτ₁₃ᶠ)
dvτdz_op = @at (Center, Center, Face) ∂z(vᶠτ₂₃ᶠ)
dwτdz_op = @at (Center, Center, Face) ∂z(wᶠτ₃₃ᶠ)

dpwdz_op = @at (Center, Center, Face) ∂z(pᶠwᶠ)

ϵ₁₁_op = @at (Center, Center, Face) 2νₑ * S₁₁_op^2
ϵ₂₂_op = @at (Center, Center, Face) 2νₑ * S₂₂_op^2
ϵ₃₃_op = @at (Center, Center, Face) 2νₑ * S₃₃_op^2
ϵ₁₂_op = @at (Center, Center, Face) νₑ * S₁₂_op^2
ϵ₁₃_op = @at (Center, Center, Face) νₑ * S₁₃_op^2
ϵ₂₃_op = @at (Center, Center, Face) νₑ * S₂₃_op^2

ω_op1 = @at (Center, Center, Center) Field(∂y(w) - ∂z(v))
ω_op2 = @at (Center, Center, Center) Field(∂z(u_E) - ∂x(w))
ω_op3 = @at (Center, Center, Center) Field(∂x(v) - ∂y(u_E))

Qᴬ_1 = Field(0.5 * (ω_op1^2))
Qᴬ_2 = Field(0.5 * (ω_op2^2))
Qᴬ_3 = Field(0.5 * (ω_op3^2))
QA = Field(Qᴬ_1 + Qᴬ_2 + Qᴬ_3)

simulation.output_writers[:vorticity_nc] =
    NetCDFOutputWriter(model, (; QA),
                       filename = "enstrophy_aligned_low_re.nc",
                       schedule = TimeInterval(output_interval))

simulation.output_writers[:velocity_nc] =
    NetCDFOutputWriter(model, (; u, v, w),
                       filename = "velocity_aligned_low_re.nc",
                       schedule = TimeInterval(output_interval))

# ---------------------------------
# Run / restart logic
# ---------------------------------

pickup_file = "./model_checkpoint_iteration108795.jld2"

if isfile(pickup_file)
    simulation.stop_time = 48.25hours
    run!(simulation, pickup = pickup_file)
else
    run!(simulation)
end
