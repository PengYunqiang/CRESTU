# CRESTU frequency-domain Rankine BEM

The solver uses the harmonic convention `exp(i*omega*t)`. Dynamic pressure
is `p=-i*omega*rho*phi`; headings are wave propagation directions measured
counterclockwise from global `+x`.

## Run

From MATLAB with `Code/Case_Wave` on the path:

```matlab
results = HydroMain();                              % Case_Wave/CRESTU.cfg
results = HydroMain('C:\path\to\another.cfg');
validation = Validate_Hemisphere_WAMIT(true);      % one-frequency real BEM
smoke = Test_MultiBody_Wave(false);                % fast 12-DOF regression
full = Test_MultiBody_Wave(true);                  % full two-body domain
symmetry = Benchmark_Symmetry(true);               % exact full/half/quarter grid
study = Run_Convergence_Study(true);               % 528/1520/3472 panels + WAMIT
```

`IPOTEN=1` assembles/factorizes the BIE and saves `<CaseName>_PotCache.mat`.
With unchanged frequencies, headings, and body-panel count, `IPOTEN=0` loads
that cache and runs force/RAO/drift post-processing without assembling the
BEM. Cache schema 3 includes symmetry parity components and the corrected
physical diffraction sign; older caches are intentionally rejected.

`PARA2` accepts either four flags (`IPOTEN IFORCE IRAD IDIFF`) or a fifth
`IDRIFT` flag. Set `IDRIFT=1` to run the direct-pressure and Kochin momentum
mean-drift modules. `PARA8 = ISX ISY` selects the retained `x>=0` and/or
`y>=0` domain. Radiation modes are grouped by reflection parity, while each
incident wave is projected into symmetric/antisymmetric components and then
recombined.

Primary outputs are `<CaseName>_Results.mat`, containing SI added mass,
damping, excitation, panel pressure, hydrostatics, and complex/amplitude/phase
RAOs. Mean drift is exported to `<CaseName>_MeanDrift.csv/.mat`. Use
`export_surface_pressure` for MAT/legacy-VTK output and
`plot_surface_pressure` for panel plots.

The free-surface outer radius is automatically enlarged so that
`R_outer-R_inner >= 1.5*lambda_max`. The radial grid has a fixed line at
`R_inner`, and the Rayleigh coefficient is selected per frequency from an
`exp(-4)` target attenuation. Production radiation damping uses the
positive-semidefinite Haskind energy-flux relation; the less robust direct
pressure value is retained as `results.diagnostics(k).pressure_damping`.

The convergence and symmetry scripts write CSV/MAT summaries and PNG plots
beside the case files. `Run_Convergence_Study(false)` and
`Benchmark_Symmetry(false)` load the latest saved reports without rerunning.

The multi-body configuration treats mass-property CG coordinates as global.
Rotational generalized normals and waterplane moments are taken about each CG.

See `UPGRADE_REPORT.md` for the validation numbers and current mean-drift
accuracy limitation.
