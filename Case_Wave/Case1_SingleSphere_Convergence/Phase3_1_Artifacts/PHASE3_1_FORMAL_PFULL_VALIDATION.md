# Phase 3.1 formal P_FULL validation

## Final gates

```text
PHASE2_3_BASELINE_REGRESSION = PASS
ABC_NORMAL_LINEAGE = PASS
PRODUCTION_ABC_UNCHANGED = YES
P_FULL_ASSEMBLY = PASS
FORMAL_FREQUENCY_SWEEP = PASS
ALL_FREQUENCIES_SOLVED = YES
MIN_RAW_RCOND = 7.316482407983828e-06
MIN_SCALED_RCOND = 1.112179795505236e-05
MAX_RELATIVE_RESIDUAL = 2.371418547490719e-14

OUTER_RADIUS_CONVERGENCE = PASS
ADOPTED_OUTER_RADIUS = 1.5R (representative-screening gate)
MAX_RADIUS_SENSITIVITY = 0.01590272065029399
MAX_FINAL_RADIUS_SENSITIVITY = 0.008271004761431481

PRODUCTION_ABC_CHANGED = NO
PRODUCTION_MESH_TOPOLOGY_CHANGED = NO
REVIEWER_BLOCKERS = 0
READY_FOR_PHASE3_2_MESH_CONVERGENCE = YES
```

## Scope and authority

The authoritative frequency list is read from `E:\CRESTU-1F_v1.0_20260824\Code\CRESTU\Case_Wave\Case1_SingleSphere_Convergence\Mesh_Fine\Case1_Fine.cfg`: 16 angular frequencies from 0.5 to 2.0 rad/s. The mesh baseline is the frozen Fine P_FULL Phase 2.3 frustum (`q_top=2.0`, `q_bottom=1.5`, body/FS radial/wall controls unchanged).

Each point was run through `run_frequency_domain_case` as an isolated clean production solve. RAW_DIAGNOSTICS was enabled only to expose the existing raw/scaled BEM condition estimates and algebraic residual records; full solved potentials were not persisted.

The current run does not modify body resolution, the FS/bottom/FARFIELD density strategy, production ABC source, WAMIT reference data, or Phase 2.3 artifacts.

## Formal sweep interpretation

All-frequency solve status: **YES**. Condition gate: **PASS**; residual gate: **PASS**; curve-continuity gate: **PASS**. The formal thresholds were rcond >= 1e-12, residual <= 1e-10, isolated curve deviation <= 0.5, and non-negligible complex phase step <= 120.0 deg.

The machine-readable response, condition/residual, and continuity records are `Phase3_1_Formal_PFull_Frequency_Sweep.csv`, `Phase3_1_Condition_Residual.csv`, and `Phase3_1_Response_Continuity_Audit.csv`.

The following coarse-grid curvature flags were resolved only by independent Phase 3.1B dense production evidence; the original coarse values are retained in the CSV:

- `A33` at omega=0.6 rad/s: dense metric 0.0245618 <= 0.5; Coarse-grid curvature resolved as dense-local continuous.
- `B33` at omega=0.6 rad/s: dense metric 0.0375286 <= 0.5; Coarse-grid curvature resolved as dense-local continuous.
- `RAO3_complex` at omega=1.4 rad/s: dense metric 0.346086 <= 0.5; Dense complex/amplitude response is continuous; finite D3 has a physical-resonance minimum.

No phase threshold, coefficient tolerance, or production response value was relaxed or smoothed.

## Full-BEM outer-radius convergence

This is a full BEM solution study, distinct from the Phase 2.3 analytic Hankel coefficient scaling. Low/mid/high, formal-response-peak, and failed continuity-audit frequencies were screened. For each frequency, the Phase 2.3 top and bottom radii were multiplied about the global origin by factors R, 1.25R, and 1.5R. The per-frequency q values were derived exactly from `q=(factor*R_baseline-r_inner)/lambda`; 2R was run only if the third level did not meet the predeclared 0.02 gate.

Changing radius necessarily adds/re-tessellates only the outer FS/bottom annuli and FARFIELD wall so that the frozen local spacing rules remain active. That derived count/hash change is the intended radius variable, not Phase 3.2 coordinated mesh convergence.

Relative sensitivities use `|Q(R2)-Q(R1)|/max(|Q(R2)|,Q_floor)` with `Q_floor=max(1e-6*formal_curve_scale,eps)`; absolute differences and the actual floor values are saved beside every comparison. `MAX_RADIUS_SENSITIVITY` covers all adjacent comparisons; the predeclared 0.02 convergence gate is applied to the final adjacent step reported as `MAX_FINAL_RADIUS_SENSITIVITY`.

## Error classification

1. **Algebraic / implementation error:** assessed by finite solves, BEM raw/scaled rcond, linear/algebraic residuals, the production one-column `D+gamma*S` oracle, source hashes, and the omega=1.5 frozen regression. Status: **PASS**.
2. **Outer-boundary truncation error:** assessed only from the full-BEM radius sensitivity table. Status: **PASS** within the representative screening set.
3. **Physical/model discrepancy:** not closed here. No WAMIT comparison was changed or used to tune the result. Any later reference discrepancy must still retain Rankine-panel integration, body/FS discretization, and source formulation as candidate causes.

## Validation and file audit

```text
IMPLEMENTED = YES
COMPILED = N/A
RUNTIME_TESTED = YES (MATLAB R2023b Update 10)
REGRESSION_PASS = PASS
PHYSICAL_VALIDATION = PARTIAL
READY_FOR_PRODUCTION = NO
```

Production dependency hashes were captured before and after the run; all 5 files were unchanged. Only the formal `Phase3_1_Artifacts` runner, CSVs, figures, checkpoint/result MAT files, source snapshots, and this report were added. Temporary per-case configs/BMFs/caches/results were created under the OS temporary directory and removed after extraction. No file was deleted from the repository.
