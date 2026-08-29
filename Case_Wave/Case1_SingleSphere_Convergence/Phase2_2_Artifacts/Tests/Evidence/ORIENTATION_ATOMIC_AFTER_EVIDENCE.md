# Double-Layer Orientation Atomic After Evidence

- Convention: `double-layer-source-normal-inward-to-fluid-v1`
- Component signs `[body, FS, bottom, outer]`: `[1 -1 1 1]`
- Orientation hash: `14308f64802f6cbb233038c51b03692c0ceb9abf16de63c6914c1fc4ed4e3948`
- Minimum stored-normal/winding alignment: `[NaN 1 1 1]`

| Test | Metric | Tolerance | Current pass | Expected after pass | Status |
|---|---:|---:|---:|---:|---|
| T2_MIXED_CANONICAL_PRODUCTION_HELPER | 5.551115123125783e-17 | 1e-10 | true | true | REPAIRED_PASS |
| T2_COMPONENT_SIGN_CONTRACT | 0 | 0 | true | true | REPAIRED_PASS |
| T2_HELPER_IN_CODE_FINGERPRINT | 0 | 0 | true | true | REPAIRED_PASS |
| T3_FS_FINITE_DEPTH_NU | 0.121590736629394 | 1e-10 | false | false | UNCHANGED_FAIL |
| T5_CYLINDER_STORED_DECAY | 0.02839456171384415 | 1e-08 | false | false | UNCHANGED_FAIL |
| T5_CYLINDER_GREEN_OUTWARD_SIGN | 1.99899337373846 | 1e-08 | false | false | UNCHANGED_FAIL |
| T5_FRUSTUM_STORED_PROJECTION | 0.1776733598773232 | 1e-08 | false | false | UNCHANGED_FAIL |
| T5_COMPLEX_KAPPA_INTERFACE | 1 | 0.5 | false | false | UNCHANGED_FAIL |
| T6_GREEN_RHS_ORIENTATION | 2.000000000000001 | 1e-12 | false | false | UNCHANGED_FAIL |
| T10_PRESSURE_TRACTION_SIGN | 1.997247499584837 | 0.02 | false | false | UNCHANGED_FAIL |
| T7_TRANSLATING_SPHERE_ADDED_MASS | 1.997247499584838 | 0.02 | false | false | UNCHANGED_FAIL |

Only T2 is repaired. Rows marked `UNCHANGED_FAIL` remain physical failures, not successful tests.
