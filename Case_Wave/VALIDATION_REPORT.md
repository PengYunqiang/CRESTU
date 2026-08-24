# CRESTU validation report

Validation date: 2026-08-24. MATLAB R2023b. All generated files are inside
`Code/Case_Wave`; WAMIT files were read only.

## Kernel and sign calibration

- Unit-square panel at `(0,0,1)`: analytic `G=0.928598`, `dG/dn=0.805432`,
  matching adaptive numerical quadrature to better than `1e-10`.
- Unit-square self source integral: `G=3.52549`, matching `4*asinh(1)`.
- Coarse 96-panel full sphere: `A11=260668 kg` versus the analytic
  `0.5*rho*volume=268344 kg`, a 2.86% error.

These checks fixed both the solid-angle kernel and the inward-fluid-normal
BIE sign (`0.5*I-D`).

## hemi_D10 versus WAMIT FullSphereIRR0

The closed Rankine domain has 300 body, 1000 free-surface, 1300 seabed, and
120 far-field panels. WAMIT `.1` values were converted from `A/rho` and
`B/(rho*omega)` to SI.

- At `omega=1 rad/s`: `A33=169199 kg` versus `155630 kg` (8.72% error);
  `B33=9120.8 kg/s` versus `90234 kg/s` (89.9% underprediction).
- At `omega=3 rad/s`: `A33=128785 kg` versus `117428 kg` (9.67% error);
  `B33=11720.9 kg/s` versus `12774.3 kg/s` (8.25% error).
- At `omega=0.05 rad/s`, Rankine `B33=-12.9 kg/s`, which is numerical zero
  at this scale. At `omega=1.5`, `B33=17759.7 kg/s`, then decreases to
  `11720.9 kg/s` at `omega=3`, confirming the sampled high-frequency decay.
- Reciprocity at `omega=1`: relative Frobenius errors are `1.29e-5` for A
  and `3.41e-5` for B.
- Hydrostatics: `Awp=78.2172 m^2`, exactly matching WAMIT's normalized
  `C33=78.2172`. The normalized roll/pitch residual is `2.829`, versus WAMIT
  `1.201`; this is cancellation-sensitive (`Ixx=486.851 m^4` and
  `volume*zb=-484.0 m^4`) and remains a mesh-convergence item.
- Excitation and RAO are finite: `|Fexc,3|=693287 N/m`, heave RAO amplitude
  `1.98763 m/m`, dynamic-system `rcond=0.129983`.

The damping mismatch at `omega=1` means this version passes sign, positivity,
reciprocity, endpoint-trend, and execution validation, but not a uniformly
converged WAMIT damping curve. Refine/tune the Rankine free-surface/radiation
boundary before treating intermediate-frequency damping as production-grade.

## Two-body 12-DOF case

- Two D=10 m hemispheres at `x=-10 m` and `x=+10 m` (20 m separation).
- Constrained two-hole free surface: 2532 panels; 600 body panels; full
  complex system `3132 x 3132`.
- Full solve completed for all 12 radiation modes and one diffraction
  heading. Added-mass/damping/hydrostatic matrices are `12 x 12`; excitation
  and RAO contain 12 rows and all values are finite.
- Assembly time was 30.435 s. RAO minimum `rcond=0.0352`.
- Multi-body symmetry errors were 1.43% (A) and 9.48% (B). This case is an
  execution/dimension benchmark, not yet a quantitative interaction benchmark.

## Cache and export

- `IPOTEN=0` reproduced cached A and excitation arrays bit-for-bit while
  skipping BEM assembly and LU.
- Cache validation covers case, frequencies, headings, environment, total
  panels, and geometry signatures.
- MAT and legacy VTK complex-pressure export passed.
