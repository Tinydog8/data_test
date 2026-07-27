# MY25 / KC04 reference profiles

Digitized eddy-viscosity profiles for validating the Mellor–Yamada 2.5 +
Kantha–Clayson (2004) Langmuir extension.

| File | Source | Columns |
|------|--------|---------|
| `kc04_fig1_KM.csv` | Kantha & Clayson (2004) Fig.1 (center panel) | `z/zi`, `KM` no-LC, `E6=1`, `E6=4` as `KM/(u★ zi)` |
| `mcwilliams1997_fig3b_KM.csv` | McWilliams et al. (1997) Fig.3b | `z/zi`, shear, Langmuir `La=0.3` as `KM/(U★ zi)` |

Hand-digitized from published figures (not author-provided tables). Absolute
pointwise error is typically a few percent of the axis range; use for
**shape / peak magnitude / E6 trend**, not sub-percent regression tests.

Compare with:

```bash
julia --project=. examples/compare_my25_kc04.jl
```

Output overlay: `output/my25_vs_kc04_fig1.csv`.
