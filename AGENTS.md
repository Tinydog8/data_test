# AGENTS.md

## Cursor Cloud specific instructions

This repository contains a single Julia module (`Mymodule.jl`) for computational fluid dynamics post-processing — specifically strain rate tensor computations on staggered grids with periodic (x/y) and solid (z) boundary conditions.

### Runtime

- **Julia 1.10+** is the only dependency. The module uses only Julia's standard library (no `Project.toml` or external packages).
- Julia is installed at `/opt/julia-1.10.4/bin/julia` (symlinked to `/usr/local/bin/julia`).

### Running / Testing

Load and test the module:

```julia
julia -e 'include("Mymodule.jl"); using .Mymodule; println("Module loaded OK")'
```

To run a full functional test with sample data:

```julia
julia -e '
include("Mymodule.jl"); using .Mymodule
u = rand(4,4,4,2); v = rand(4,4,4,2); w = rand(4,4,4,2)
∂x(u, 1.0); ∂y(u, 1.0); ∂z(u, 1.0)
InterΣᶜᶜᶜ(u,v,w,1.0,1.0,1.0)
ΣΣ(u,v,w,1.0,1.0,1.0)
trΣ(u,v,w,1.0,1.0,1.0)
Σᵢ₃ᶜᶜᶜ(u,v,w,1.0,1.0,1.0)
println("All functions OK")
'
```

### Notes

- There are no lint tools, test frameworks, or build steps configured in this repo.
- The module has no external dependencies — `include("Mymodule.jl")` is all that's needed.
- Comments in the source are in Chinese and describe staggered-grid interpolation strategies for avoiding boundary issues.
