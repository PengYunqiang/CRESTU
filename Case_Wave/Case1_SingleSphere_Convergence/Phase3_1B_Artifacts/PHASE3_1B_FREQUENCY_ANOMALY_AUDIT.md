# Phase 3.1B frequency-response anomaly root-cause audit

## Final decision

The two Phase 3.1 blockers are closed without changing the production ABC,
Rankine formulation, mesh-topology generator, A/B extraction formula, or
production response values.

```text
PHASE2_3_BASELINE_REGRESSION = PASS
PRODUCTION_ABC_CHANGED = NO
PRODUCTION_MESH_TOPOLOGY_CHANGED = NO

OMEGA_0P6_FIRST_ANOMALOUS_LAYER = POTENTIAL
RAW_RADIATION_FORCE_CONTINUITY = PASS
AB_EXTRACTION_CONTINUITY = PASS

ANOMALY_PRESENT_WITH_REMESH = NO
ANOMALY_PRESENT_WITH_FIXED_GEOMETRY = NO

FREQUENCY_LOCAL_MESH_DISCONTINUITY = NO
FREQUENCY_LOCAL_DISCRETIZATION_CAUSAL = NO

RANKINE_QUADRATURE_SENSITIVITY = UNRESOLVED
RADIATION_BOUNDARY_SENSITIVITY = NO
TRUNCATED_DOMAIN_NUMERICAL_MODE = UNRESOLVED

OMEGA_0P6_ROOT_CAUSE = coarse frequency-grid curvature aliasing of a dense-local continuous Rankine-only response

OMEGA_1P4_RAO3_CLASSIFICATION = PHYSICAL_RESONANCE

PRODUCTION_FIX_IMPLEMENTED = NO
PRODUCTION_FIX_DESCRIPTION = no production numerical formula or mesh topology changed
VALIDATION_HARNESS_FIX_IMPLEMENTED = YES (dense-evidence adjudication only)

FORMAL_FREQUENCY_SWEEP_AFTER_FIX = PASS
OUTER_RADIUS_REGRESSION = PASS
FREQUENCY_RESPONSE_CONTINUITY = PASS

REVIEWER_BLOCKERS = 0
READY_FOR_PHASE3_2_MESH_CONVERGENCE = YES
```

`OMEGA_0P6_FIRST_ANOMALOUS_LAYER=POTENTIAL` means that the wide/coarse
curvature signal is already present in the solved body-potential vector; it
does not mean that the dense potential is discontinuous. On the required
0.59/0.60/0.61 triplet, every inspected layer is continuous.

## Production radiation lineage and convention

The executed production path uses the harmonic convention
`exp(+i*omega*t)` and a unit generalized displacement. For heave:

```text
q3 = i*omega*n3
rhs = -S*q3
p = -i*rho*omega*phi3
F33_rad = i*omega*rho*sum(n3*phi3*dS)
A33 = real(F33_rad)/omega^2
B33 = -imag(F33_rad)/omega
F33_rad = omega^2*A33 - i*omega*B33
```

The diagonal reciprocal projection cannot change A33 or B33. The reconstructed
panel-force sum matches `results.diagnostics.radiation_force(3,3)` with a
maximum relative difference of `4.74e-16`. No frequency-dependent branch was
found in the A/B extraction, density scaling, or omega normalization.

The complete 0.58--0.62 body and full-solution phi3 payload is stored in
`Phase3_1B_Radiation_Phi3_Lineage.mat`; its scalar probes, norms, raw force,
and coefficient lineage are in `Phase3_1B_Radiation_Quantity_Lineage.csv`.

## Q1: omega=0.6 A33/B33

### Layer localization

The wide/coarse curvature is visible before force post-processing, but it
vanishes under the required dense-local test:

| Layer/quantity | Wide or formal deviation | Dense 0.59/0.60/0.61 deviation | Gate |
|---|---:|---:|---:|
| body phi3 complex vector | 0.638471 (0.54/0.60/0.66) | 0.036302 | PASS |
| raw F33 complex force | 0.583604 (0.54/0.60/0.66) | 0.032938 | PASS |
| A33 | 0.734631 (formal 0.5/0.6/0.7) | 0.024562 | PASS |
| B33 | 0.596420 (formal 0.5/0.6/0.7) | 0.037529 | PASS |

The apparent isolated dip/spike therefore came from applying a three-point
linear-curvature gate to a response segment sampled only every 0.1 rad/s. It
is not an extraction bug and not a dense-local solver discontinuity.

### Frequency-local discretization

At 0.59, 0.60, and 0.61 rad/s, the actual production mesh has identical
discrete controls:

```text
N_body/N_FS/N_bottom/N_FARFIELD/N_total = 588/680/1164/576/3008
legacy/top/bottom/common theta counts = 56/72/56/72
top/bottom radial layers = 2/1
wall vertical layers = 8
```

Top/bottom radii and spacings vary monotonically with wavelength, so the exact
component and combined geometry hashes correctly differ at every frequency.
That expected coordinate/hash change is not a discrete jump. The first nearby
integer change occurs at 0.64 rad/s, where the wall changes from eight to seven
layers; it is outside the 0.59/0.60/0.61 target triplet.

The omega=0.60 fixed-geometry run reuses one exact combined geometry hash at
0.56, 0.58, 0.59, 0.60, 0.61, 0.62, and 0.64 rad/s. Its local raw-F/A/B
deviations at 0.60 are 0.1394/0.2045/0.1120, all below 0.5. Thus:

```text
FREQUENCY_LOCAL_DISCRETIZATION_CAUSAL = NO
```

The fixed and remeshed coefficient magnitudes are not identical away from the
anchor; frequency-local outer geometry still has quantitative discretization
sensitivity that belongs in Phase 3.2. It is not the cause of an isolated
0.60 discontinuity.

### Lower-priority candidates

- The R/1.25R/1.5R 0.60 screen passes the 2% gate. The R->1.25R relative
  changes are 0.00648 (A33) and 0.01590 (B33); the 1.25R->1.5R changes are
  0.00626 and 0.000675. Radiation-boundary sensitivity is therefore `NO`
  within this declared screening scope.
- Rankine-panel quadrature perturbation was not triggered because neither
  remeshed nor fixed geometry retained a dense-local anomaly. The production
  Rankine panel integral has no frequency-dependent quadrature-order switch,
  so this status remains `UNRESOLVED`, not a claimed cause.
- A truncated-domain numerical-mode experiment was not triggered for the same
  reason. Rcond/residual trends are smooth and the radius screen passes, but an
  explicit moving-mode/singular-vector study was unnecessary; status is
  `UNRESOLVED` rather than a positive attribution.
- No lid, extended BIE, irregular-frequency suppression, or source-formulation
  change was implemented. No evidence supports the classical irregular-
  frequency mechanism of free-surface Green-function BEM.

### Read-only reference diagnostic

The existing WAMIT IRR3 reference is smooth around 0.6. At 0.60 rad/s it gives
A33=212307 kg and B33=40196 kg/s, whereas the current Rankine result gives
66111 kg and 73275 kg/s. This supports a Rankine-only broad model/discretization
feature, not an isolated WAMIT-like spike. The off-grid 0.58 and 0.62 reference
values in the comparison CSV are explicitly labeled PCHIP diagnostic
interpolations and were not used to tune any parameter.

## Q2: omega=1.4 RAO3

The production heave dynamic matrix is

```text
D = -omega^2*(M+A) + i*omega*B + C
RAO = D\F_excitation
```

The 1.30:0.02:1.50 sweep shows a continuous finite resonance:

```text
omega at minimum |D3| = 1.40 rad/s
D3(1.40) = 13970.45 + 70701.28i N/m
|D3(1.40)| = 72068.33 N/m
RAO3(1.40) = 1.55798 - 2.34879i
max |RAO3| = 2.81853 at 1.40 rad/s
max dense normalized complex step = 0.34609
max dense amplitude relative step = 0.26132
phase amplitude protection triggered = NO
```

Real(D3) changes sign between 1.40 and 1.42 while damping keeps the denominator
finite. Amplitude, complex response, and unwrapped phase rotate continuously.
WAMIT also shows a smooth resonance in this band. The classification is
`PHYSICAL_RESONANCE`, not a near-zero phase artifact or numerical
discontinuity.

The old failure was not a raw phase-step failure: its phase step was only
47.79 degrees against a 120-degree limit. It was the coarse 1.3/1.4/1.5
three-point complex-curvature metric. The harness now retains that preliminary
failure and requires the dense complex/amplitude/D3 evidence before changing
only the validation status.

## Minimal changes

1. `run_frequency_domain_case` gained default-off diagnostic name-value
   controls for a physics-frequency override, an exact fixed-geometry anchor,
   and skipping unrelated off-surface diagnostics. Ordinary one-argument
   production behavior and all numerical formulas are unchanged.
2. The Phase 3.1 validation runner now persists actual radius counts,
   spacings, and component/combined hashes; writes a dedicated frozen-
   regression difference CSV; and binds formal/radius checkpoints to source,
   config, geometry-control, frequency-list, validation-option, runner-version,
   and runner-hash fingerprints.
3. Coarse continuity failures may be adjudicated only by the required dense
   evidence. Preliminary failures and all production values remain visible.
4. `Verify_Phase3_1B_Final_Evidence.m` is the authoritative final gate. It
   verifies both runner hashes, the current 41-file production fingerprint,
   config hashes, exact frequency vectors, evidence-file hashes, target-triplet
   discrete controls, conditioning/residual limits, radius sensitivities, and
   the independent reviewer gate. It supersedes the historical
   `Phase3_1B_Provisional_Decisions.csv`, whose unique-hash heuristic incorrectly
   labeled continuous frequency-local geometry changes as a discontinuity.

No coefficient smoothing, interpolation replacement, frequency special case,
ABC change, mesh-topology change, WAMIT modification, or physical-formulation
change was made.

## Validation actually run

- MATLAB R2023b Update 10, real installed runtime.
- Remeshed local sweep: 0.54, 0.56, 0.58, 0.59, 0.60, 0.61, 0.62, 0.64,
  0.66 rad/s.
- Exact 0.60 fixed-geometry sweep: 0.56, 0.58, 0.59, 0.60, 0.61, 0.62,
  0.64 rad/s.
- RAO/D3 sweep: 1.30:0.02:1.50 rad/s.
- Clean formal 16-frequency sweep: 0.5:0.1:2.0 rad/s.
- Frozen Phase 2.3 omega=1.5 regression: complex F3 relative difference
  4.28e-13; raw/scaled rcond and panel counts exactly reproduced.
- Outer-radius regression at 0.5, 0.6, 1.2, 1.4, and 2.0 rad/s, each at
  R/1.25R/1.5R: PASS.
- Formal minimum raw/scaled rcond: 7.316e-6 / 1.112e-5; maximum linear or
  algebraic relative residual: 2.371e-14.
- Fingerprint-matched checkpoint reuse smoke test: PASS with zero BEM reruns.
- Final evidence verifier: PASS.
- Independent read-only reviewer: 12 scopes, zero blockers.

## Validation and file-hygiene audit

```text
IMPLEMENTED = YES
COMPILED = N/A
RUNTIME_TESTED = YES (MATLAB R2023b Update 10)
REGRESSION_PASS = YES
PHYSICAL_VALIDATION = PARTIAL
REVIEWER_BLOCKERS = 0
READY_FOR_PRODUCTION = NO
READY_FOR_PHASE3_2_MESH_CONVERGENCE = YES
```

`PHYSICAL_VALIDATION=PARTIAL` is deliberate: the current task isolates the two
frequency-response blockers and uses WAMIT only as a diagnostic oracle. It does
not complete Phase 3.2 mesh convergence or full physical validation.

No repository file was deleted. No temporary/debug directory was added to the
project. Runtime configs, BMFs, caches, raw diagnostics, and result MAT files
were created only in the OS temporary directory and cleaned after extraction.
The pre-existing dirty worktree was preserved; unrelated user changes were not
reverted or overwritten.
