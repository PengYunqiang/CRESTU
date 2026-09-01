# Phase 2.3 P_FULL ABC normal/operator audit

## Final status

```text
ABC_NORMAL_LINEAGE = PASS
ABC_WALL_NTHETA_STATUS = RESOLVED
ABC_IMPLEMENTATION_ERROR = PASS
HANKEL_ASYMPTOTIC_ERROR_CLASSIFIED = YES
P_FULL_ASSEMBLY = PASS
OMEGA_1P5_FOCUSED_CASE = PASS
FORMAL_FREQUENCY_SWEEP = NOT_RUN
READY_FOR_FORMAL_PFULL_FREQUENCY_SWEEP = YES
REVIEWER_BLOCKERS = 0
```

Decision rule **A** applies. The active production FARFIELD wall has
`max|n_theta| < 1e-10`; no azimuthal, modal, tangential, or nonlocal ABC
operator is required. No production ABC or mesh-topology source was changed
in Phase 2.3.

## 1. Scope and runtime authority

The audit used MATLAB R2023b Update 10 and the focused full-domain case
`Frustum_Top_Only.cfg` at `omega=1.5 rad/s`, `q_top=2`, `q_bottom=1.5`.
The current production geometry was rebuilt in memory through the same mesh
generators and `apply_rankine_outer_truncation` used by the frequency-local
production path. The audit did not use a cached matrix or reload a generated
FARFIELD BMF for assembly.

The runtime path audit resolved all five critical functions to the formal
`Source Code` tree. None resolved to `_Phase2_2_Pre*` shadow directories; see
`Phase2_3_Runtime_Function_Resolution.csv`.

## 2. Normal lineage and the source of `0.016359`

### Active production lineage

The production lineage is:

```text
build_frequency_local_domain / focused in-memory equivalent
  -> apply_rankine_outer_truncation
  -> FARFIELD mesh.normals (diagonal-cross winding normal)
  -> merge_domain_geometry copies normals unchanged
  -> get_rankine_outer_absorbing_bc selects tail FARFIELD rows by component counts
  -> assemble_rankine_matrix consumes the same merged normals and gamma
```

For this case the merged component counts are body `588`, free surface `680`,
bottom `1212`, and FARFIELD `272`, so the production ABC consumes global
source panels `2481:2752`. The positional tail mask exactly equals
`find(panel_type==5)`. The FS/bottom count-transition annuli remain boundary
types `2` and `4` and are not selected.

`merge_domain_geometry` does not recalculate or project the FARFIELD normals.
The cross-product normal is used by the wall generator and copied into the
merged normal array. The analytic mother-line normal is used only as an
independent geometry oracle.

### Four normal definitions on the active wall

All cylindrical bases use the requested panel centroid definition
`theta_c=atan2(y_c,x_c)`.

| Normal definition | `max|n_theta|` | Role |
|---|---:|---|
| production merged/stored normal | `1.99840144432528e-15` | actual assembly input |
| diagonal cross-product geometric normal | `1.99840144432528e-15` | independent reconstruction |
| analytic axisymmetric mother-line normal | `1.11022302462516e-16` | continuous geometry oracle |
| legacy benchmark first-triple algorithm applied to production panels | `2.60902410786912e-15` | definition comparison |

The merged/component normal difference is exactly zero. Minimum stored/cross
alignment is `0.99999999999999967`; minimum stored/benchmark-algorithm
alignment is `0.99999999999999978`. The structured seam audit covers every
layer, the final-to-first circumferential connection, and adjacent-layer
connectivity; it passed.

`get_rankine_outer_absorbing_bc` uses the algebraically equivalent Cartesian
centroid projection and reports `max|n_theta|=1.94719585036681e-15` during
the full assembly. Its maximum arithmetic difference from the explicit
`atan2` basis projection is `4.17393574781679e-16`.

The requested 20-panel active-wall detail is in
`Phase2_3_Normal_Lineage_Top20.csv`. The node IDs in that CSV are
coordinate-tolerance equivalence IDs reconstructed for connectivity audit;
BMF does not carry native node IDs.

### Final origin of `0.016359`

The Phase 2.2 Hankel script did not load the current production wall. It
hard-coded the frozen file `Frustum_Fine_Top_Only_farfield.bmf`, whose header
is `Nz=5 Ntheta=68->56 unequal-zipper` and whose panel count is `484`. It also
hard-asserted 484 panels and hard-coded `0.01635907986834724` as a synthetic
test magnitude.

Re-running its unique-vertex-centroid and first-noncollinear-triple normal
algorithm gives `max|n_theta|=0.016359080607805082` on that stale file. The
maximum is tied, to serialization precision, at old local panels
`373, 401, 429, 457` (old merged global sources `2789, 2817, 2845, 2873`).
They are three-unique-node zipper triangles in the old `60 -> 56` transition
layer. The top 20 legacy panels are recorded in
`Phase2_3_Legacy_Hankel_0p016359_Top20.csv`.

Therefore:

- it was not a cache or shadow result;
- it was not caused by node theta (the old benchmark also used panel-centroid
  theta);
- it did not select the current FS/bottom transition annuli;
- it did select an obsolete FARFIELD wall whose count transitions were on
  the wall itself;
- it is a faceted unequal-zipper artifact, not the normal of the current
  axisymmetric target wall.

The two quoted old values differ by about `7.39e-10` because the in-memory old
value was written to a nine-digit BMF and then reconstructed.

The Phase 2.2 equal-count delivery BMF is also a nine-digit text snapshot;
reloading that snapshot gives `max|n_theta|=1.23267823726181e-9`. It is
explicitly classified as `DELIVERY_SNAPSHOT_NOT_RELOADED`, not as the active
assembly input. Production rebuilds and assembles the in-memory wall, whose
authoritative value is `1.99840144432528e-15`.

## 3. Hankel operator error decomposition

The production continued complex wavenumber `kappa` is used for the active
benchmark. For mode `m`,

```text
Phi = H_m^(2)(kappa*r) exp(i*m*theta) cosh(kappa*(z+h))/cosh(kappa*h)
q_full = n_r*d_r Phi + (n_theta/r)*d_theta Phi + n_z*d_z Phi
q_mer  = n_r*d_r Phi + n_z*d_z Phi
q_abc  = [(-n_r)*(i*kappa+1/(2r))
          +n_z*kappa*tanh(kappa*(z+h))] Phi
q_actual = outerBC.gamma .* Phi
```

The three errors are:

```text
E_THETA          = q_full - q_mer
E_ASYMPTOTIC     = q_mer - q_abc
E_IMPLEMENTATION = q_actual - q_abc
```

Pointwise relative metrics use a field-scale denominator based on both
references and `|kappa*Phi|`, with a dimensionally consistent global floor.
The L2 metric is panel-area weighted. This avoids false amplification where
the exact projected derivative is close to zero.

| m | Error | Max relative | P95 relative | Weighted L2 relative |
|---:|---|---:|---:|---:|
| 0 | `E_THETA` | `0` | `0` | `0` |
| 0 | `E_ASYMPTOTIC` | `2.06141618710528e-4` | `2.06141618710395e-4` | `1.57251389827245e-4` |
| 0 | `E_IMPLEMENTATION` | `4.45439720188121e-16` | `2.38680489195246e-16` | `1.34881104383973e-16` |
| 1 | `E_THETA` | `1.31096241888527e-16` | `6.55481209442632e-17` | `2.49420801665328e-17` |
| 1 | `E_ASYMPTOTIC` | `6.18142745654004e-4` | `6.18142745653815e-4` | `4.71598545584704e-4` |
| 1 | `E_IMPLEMENTATION` | `3.79350815437889e-16` | `2.70340772893989e-16` | `1.41739300345091e-16` |

Across all radius factors and both modes, coefficient/formula
`E_IMPLEMENTATION` has maximum relative error `4.65028919029354e-16` and
maximum weighted-L2 relative error `1.48859988789785e-16`. Detailed absolute
and relative metrics are in `Phase2_3_Hankel_Error_Decomposition.csv`.

## 4. Radius scaling and error classification

Radius scaling is analytic only: the production mesh and ABC are not
modified. Panel theta, z, normal, and quadrature weights are held fixed while
the cylindrical radius is evaluated at `R, 1.25R, 1.5R, 2R`. A canonical
inward cylinder coefficient is reported separately so the fitted order is not
contaminated by the sloped-wall vertical projection.

For the production complex `kappa`, the canonical fitted decay orders are
`1.966697` for `m=0` and `1.966410` for `m=1`. The active sloped-wall
weighted-L2 orders are `1.939579` and `1.939201`. Both are consistent with the
expected higher-order residual close to `R^-2`.

The quoted legacy `0.17843%` value is reproduced with the old wall and its
real physical `k=0.229436147969854 1/m`. Its finite-radius `m=1`
`E_ASYMPTOTIC` scales as:

| Radius | Max relative error | Percent |
|---|---:|---:|
| `R` | `1.78430557460698e-3` | `0.1784305575%` |
| `1.25R` | `1.14356698626774e-3` | `0.1143566986%` |
| `1.5R` | `7.94757843018731e-4` | `0.0794757843%` |
| `2R` | `4.47397284555022e-4` | `0.0447397285%` |

The fitted old-wall physical-k orders are `1.989963` for `m=0` and `1.991279`
for `m=1`. This is asymptotic truncation/model error, not software
implementation error. Full results are in
`Phase2_3_Hankel_Radius_Scaling.csv` and
`Phase2_3_Legacy_0p17843_Radius_Scaling.csv`.

## 5. Production assembly and focused case

Only after the four pre-assembly gates and two independent reviewer scopes
reported zero blockers was the full assembly started.

The production assembly completed for `N=2752`, with `272` FARFIELD source
columns. The implementation oracle independently reconstructed every
FARFIELD column as

```text
A_FF_oracle = D_canonical_FF + S_FF .* gamma_formula'
```

using the same full-domain parity, production source-orientation signs,
`-1/(4*pi)` double-layer convention, `+0.5` self term, single-layer scaling,
and an independently evaluated analytic gamma from the actual centers,
normals, and production `kappa`.

| Assembly implementation metric | Value |
|---|---:|
| gamma maximum relative error | `1.60306799775579e-16` |
| 272-column maximum stable residual | `2.00150994626751e-16` |
| 272-column P95 stable residual | `1.42871431869814e-16` |
| full FARFIELD block Frobenius residual | `5.19521468646348e-17` |
| production assembly time | `46.994819 s` |
| all-column oracle time | `2.418316 s` |

The per-column evidence is in
`Phase2_3_PFULL_All_Farfield_Column_Oracle.csv`.

This oracle shares the already-used `rankine_panel_integrals` kernel with the
production assembly. It independently verifies FARFIELD selection, source
orientation, signs, self terms, gamma wiring, and all-column assembly; it is
not an independent validation of the Green-function quadrature kernel, which
is outside this phase's scope.

The single `omega=1.5 rad/s` focused diffraction solve passed with linear and
algebraic relative residual `8.00262089350084e-14`. Raw and scaled reciprocal
condition estimates were `2.31950141753968e-5` and
`6.18624988260384e-5`. This is runtime and algebraic evidence, not a WAMIT or
experimental physical validation.

## 6. Decision

```text
PRODUCTION_ABC_MODIFIED = NO
NEW_DTHETA_OPERATOR_REQUIRED = NO
DECISION_RULE = A
FORMAL_FREQUENCY_SWEEP_EXECUTED = NO
READY_FOR_FORMAL_PFULL_FREQUENCY_SWEEP = YES
```

The current local first-order production ABC is implemented consistently
with its own analytic formula to machine precision. Its exact-Hankel
finite-radius discrepancy is correctly classified as first-order asymptotic
model truncation. It must not be relabeled as an implementation failure.

## 7. Validation and file audits

```text
IMPLEMENTED = YES
COMPILED = N/A
MATLAB_R2023B_RUNTIME_TESTED = YES
CHECKCODE_PASS = YES (0 findings on both Phase 2.3 runners)
REGRESSION_PASS = YES (Phase 2.3 scoped analytic/assembly/focused gates)
NUMERICALLY_VERIFIED = YES (within the listed analytic and matrix-oracle scope)
PHYSICAL_VALIDATION = PARTIAL (analytic Hankel only; no WAMIT comparison)
REVIEWER_BLOCKERS = 0
READY_FOR_PRODUCTION = NO
READY_FOR_FORMAL_PFULL_FREQUENCY_SWEEP = YES
```

File hygiene audit:

- Added directory: `Phase2_3_Artifacts` only.
- Production source files modified by Phase 2.3: none.
- Phase 2.2 artifacts overwritten: none.
- Deleted files: none.
- Temporary MATLAB preference directory: removed after validation.
- Formal additions: two reproducible MATLAB runners, CSV audit evidence, two
  compact MAT result records, reviewer gate, and this report. Every added
  file is Phase 2.3 evidence intended to remain with the phase delivery.

The repository contained pre-existing Phase 2.2 modified/untracked work when
Phase 2.3 began. It was preserved and not reformatted, reverted, or moved.
