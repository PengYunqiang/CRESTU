# Phase 2.2 Analytic Fail-Before Evidence

This is a physical fail-before audit. `currentPass=false` remains a failure; an expected failure is not counted as a validated test.

- Time convention: `exp(+i*omega*t)`
- Fine body: `E:\CRESTU-1F_v1.0_20260824\Code\CRESTU\Case_Wave\Case1_SingleSphere_Convergence\Mesh_Fine\hemi_D10_fine_body.bmf` (588 panels)
- CSV: `E:\CRESTU-1F_v1.0_20260824\Code\CRESTU\Case_Wave\Case1_SingleSphere_Convergence\Phase2_2_Artifacts\Tests\Evidence\Analytic_FailBefore_Results.csv`
- MAT: `E:\CRESTU-1F_v1.0_20260824\Code\CRESTU\Case_Wave\Case1_SingleSphere_Convergence\Phase2_2_Artifacts\Tests\Evidence\Analytic_FailBefore_Results.mat`

| Test | Category | Before expectation | Current pass | Relative error | Tolerance | Gate |
|---|---|---:|---:|---:|---:|---|
| T1_TIME_PHASE | time/amplitude | EXPECTED_PASS_BEFORE | true | 1.46031421e-16 | 1e-12 | ACTIVE |
| T1_PEAK_AMPLITUDE_LINEARITY | time/amplitude | EXPECTED_PASS_BEFORE | true | 0 | 1e-12 | ACTIVE |
| T2_INWARD_CLOSED_CONTROL | Green normal/winding | EXPECTED_PASS_BEFORE | true | 5.55111512e-17 | 1e-10 | ACTIVE |
| T2_INWARD_WINDING_CONTROL | Green normal/winding | EXPECTED_PASS_BEFORE | true | 0 | 1e-12 | ACTIVE |
| T2_CURRENT_MIXED_ORIENTATION | Green normal/winding | EXPECTED_FAIL_BEFORE | false | 0.217952892 | 1e-10 | ACTIVE |
| T2_MIXED_WINDING_IS_SELF_CONSISTENT | Green normal/winding | EXPECTED_PASS_BEFORE | true | 0 | 1e-12 | ACTIVE |
| T2_FINE_BODY_WINDING_ALIGNMENT | Green normal/winding | EXPECTED_PASS_BEFORE | true | 3.33066907e-16 | 1e-12 | ACTIVE |
| T3_FS_FINITE_DEPTH_NU | free-surface Robin | EXPECTED_FAIL_BEFORE | false | 0.121590737 | 1e-10 | ACTIVE |
| T4_BOTTOM_NEUMANN | bottom Neumann | EXPECTED_PASS_BEFORE | true | 0 | 1e-12 | ACTIVE |
| T5_CYLINDER_STORED_DECAY | outer outgoing | EXPECTED_FAIL_BEFORE | false | 0.0283945617 | 1e-08 | ACTIVE |
| T5_CYLINDER_GREEN_OUTWARD_SIGN | outer outgoing | EXPECTED_FAIL_BEFORE | false | 1.99899337 | 1e-08 | ACTIVE |
| T5_FRUSTUM_STORED_PROJECTION | outer outgoing | EXPECTED_FAIL_BEFORE | false | 0.17767336 | 1e-08 | ACTIVE |
| T5_COMPLEX_KAPPA_INTERFACE | outer outgoing | EXPECTED_FAIL_BEFORE | false | 1 | 0.5 | NOT_IMPLEMENTED |
| T6_BODY_NEUMANN_VALUE | radiation RHS | EXPECTED_PASS_BEFORE | true | 0 | 1e-12 | ACTIVE |
| T6_GREEN_RHS_ORIENTATION | radiation RHS | EXPECTED_FAIL_BEFORE | false | 2 | 1e-12 | ACTIVE |
| T10_SPHERE_QUADRATURE_CONTROL | pressure traction | EXPECTED_PASS_BEFORE | true | 0.00275250042 | 0.02 | ACTIVE |
| T10_PRESSURE_TRACTION_SIGN | pressure traction | EXPECTED_FAIL_BEFORE | false | 1.9972475 | 0.02 | ACTIVE |
| T7_TRANSLATING_SPHERE_ADDED_MASS | pressure A/B extractor | EXPECTED_FAIL_BEFORE | false | 1.9972475 | 0.02 | ACTIVE |
| T7_TRANSLATING_SPHERE_DAMPING | pressure A/B extractor | EXPECTED_PASS_BEFORE | true | 0 | 1e-12 | ACTIVE |

## Interpretation

- Passing control rows establish that the analytic oracles and existing Fine sphere mesh are active.
- Failing rows identify convention defects before any production edit.
- The complex-kappa row is `NOT_IMPLEMENTED`; it is not evidence of a validated complex radiation root.
