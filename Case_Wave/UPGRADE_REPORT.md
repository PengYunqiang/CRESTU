# CRESTU symmetry, mean-drift, and convergence upgrade

Validation date: 2026-08-24, MATLAB R2023b. All generated/modified files are
inside `Code/`. The supplied manuals and `WAMIT/FullSphereIRR0`/`IRR3` were
read only.

## Implemented physics and solver changes

- `ISX`/`ISY` geometry reduction retains `x>=0` and/or `y>=0`. Rankine source
  and doublet kernels superpose 2 or 4 reflected panels. Polar-vector and
  axial-vector parities are applied to all `6N` modes.
- Incident waves are projected into each active symmetry parity using their
  reflected headings. The parity solutions are recombined for physical
  diffraction pressure and generalized forces.
- Cache schema 3 records parity components and rejects caches with the old
  diffraction sign convention.
- The diffraction algebraic RHS sign was corrected. At `omega=1`, the fine
  direct heave excitation changed from 682.7 kN/m to 411.1 kN/m; the latter
  agrees with the independent Haskind result and WAMIT (418.0 kN/m).
- The free-surface grid now puts row `NR_near+1` exactly at `R_inner` instead
  of a hard-coded 55% radius. `R_outer-R_inner` is expanded to at least
  `1.5*lambda_max`; a frequency-dependent quadratic Rayleigh profile targets
  `exp(-4)` attenuation.
- Radiation damping used by RAO is evaluated from the Haskind energy-flux
  identity and is symmetric positive semidefinite. Direct pressure damping is
  retained in diagnostics. This removes the cancellation-sensitive
  mid-frequency underprediction without empirical WAMIT scaling.
- `6.MeanDriftLoads/` implements the four requested direct-pressure terms,
  local least-squares velocity/Hessian recovery, Kochin momentum loads, full
  symmetry-field reconstruction, WAMIT `.8/.9` ingestion, and CSV/MAT export.

## Exact-grid symmetry benchmark

The half and full reference grids are reflections of the same quarter grid,
so all cases represent exactly 740 full-domain panels.

| Domain | Matrix unknowns | Runtime (s) | Speedup | max relative result difference |
|---|---:|---:|---:|---:|
| Full | 740 | 2.413 | 1.000 | 0 |
| Half-X | 370 | 2.283 | 1.057 | 2.3e-16 |
| Quarter-XY | 185 | 2.234 | 1.080 | 8.3e-15 |

The unknown count falls by 2x/4x. Overall speedup is modest for this small
case because reflected kernel evaluation dominates assembly; dense LU cost
scales more favorably at larger matrices. `A33`, energy `B33`, and complex
heave excitation agree to machine precision across the three domains.

## Grid convergence against WAMIT IRR0

The full-domain grids contain 528, 1520, and 3472 total boundary panels.
Fine-grid relative errors are:

| omega (rad/s) | A33 | B33 | heave excitation | heave RAO |
|---:|---:|---:|---:|---:|
| 0.75 | 23.27% | 10.19% | 4.86% | 0.30% |
| 1.00 | 1.65% | 2.71% | 1.65% | 1.24% |
| 1.25 | 0.50% | 7.82% | 5.50% | 4.02% |

At the requested anomaly frequency, fine `B33=87,785 kg/s` versus WAMIT
`90,234 kg/s`. The dedicated single-frequency regression with a 1.5-wave
sponge gives `94,637 kg/s` (4.88% error). IRR0 and IRR3 reference curves are
both included in `Convergence_Study.png`.

The coarse grid is intentionally under-resolved and produces a nonphysical
negative `A33` at `omega=1`; it is retained to make the convergence failure
visible rather than silently filtering it.

## Mean-drift validation status

Fine-grid surge coefficients (`Cd=F/(0.5*rho*g*A^2*L)`) are:

| omega | CRESTU near | CRESTU far | WAMIT pressure (.9) | WAMIT momentum (.8) |
|---:|---:|---:|---:|---:|
| 0.75 | 0.05787 | 0.05960 | 0.0000063 | -0.0000265 |
| 1.00 | 0.35428 | 0.31954 | 0.01642 | 0.01398 |
| 1.25 | 2.51496 | 1.29471 | 0.63345 | 0.63827 |

Near/far methods have the same trend and are within 10% at `omega=1`, but
both still overpredict WAMIT materially. The module is execution-complete and
auditable term-by-term, but is **not yet quantitatively validated for
production drift loads**. Likely remaining sensitivities are the low-order
surface Hessian and the Kochin normalization/cross term. The plot deliberately
shows this gap (`MeanDrift_Validation.png`).

## Modified code

- `1.Input/read_config.m`
- `2.Mesh/build_bmf_domain.m`, `extract_waterline.m`,
  `generate_free_surface_bmf.m`
- `4.Potential/assemble_rankine_matrix.m`, `load_potential_cache.m`,
  `solve_radiation_freq.m`
- `5.Force/compute_hydrodynamic_coeffs.m`,
  `compute_hydrostatic_matrix.m`
- `Case_Wave/HydroMain.m`, `run_frequency_domain_case.m`, `README.md`

## New code and cases

- Symmetry/mesh: `reduce_mesh_by_symmetry.m`, `expand_mesh_by_symmetry.m`,
  `complete_waterline_by_symmetry.m`, `generate_reduced_seabed_mesh.m`,
  `tune_sponge_layer.m`
- Symmetry/potential: `get_mode_parities.m`, `symmetry_force_weights.m`,
  `decompose_incident_wave_symmetry.m`, `solve_wave_dispersion.m`
- Force/WAMIT: `compute_haskind_excitation.m`,
  `compute_radiation_damping_energy.m`, `read_wamit_excitation.m`,
  `read_wamit_rao.m`
- Mean drift: every file in `6.MeanDriftLoads/`
- Drivers: `Benchmark_Symmetry.m`, `Run_Convergence_Study.m`,
  `CRESTU_Sym*.cfg`, `CRESTU_Convergence_*.cfg`, and
  `CRESTU_DampingFix*.cfg`
- Generated validation artifacts: `Symmetry_Benchmark.csv/.mat`,
  `Convergence_Study.csv/.mat/.png`, `MeanDrift_Validation.png`, plus the
  per-grid potential/result/drift files.

## Run commands

```matlab
cd('F:\...\Rankine源\Code')
addpath(genpath(pwd))

symmetry = Benchmark_Symmetry(true);
study = Run_Convergence_Study(true);
results = HydroMain(fullfile('Case_Wave','CRESTU_Convergence_Fine.cfg'));

% Read saved reports without solving again:
symmetry = Benchmark_Symmetry(false);
study = Run_Convergence_Study(false);
```

Set the fifth `PARA2` flag to `1` for drift. For cache-only post-processing,
set `IPOTEN=0` without changing case name, geometry, frequencies, headings,
or symmetry flags.
