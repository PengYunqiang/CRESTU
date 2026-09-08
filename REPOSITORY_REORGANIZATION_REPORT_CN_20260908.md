# CRESTU Repository Reorganization Report

Date: 2026-09-08  
Release branch during reorganization: `codex/release-20260908`  
NextStep policy: `DO_NOT_TOUCH_NEXTSTEP = TRUE`

## 1. Result

The repository is now a small, portable MATLAB release tree rather than the full research working directory.

```text
HANDOFF_READY = YES
MAIN_WORKTREE_CLEAN = YES
FRESH_CLONE_TEST = PASS
MATLAB_SETUP_TEST = PASS
MINIMAL_CASE_TEST = PASS
VALIDATION_SMOKE_TEST = PASS
```

The minimal case was executed in real MATLAB R2023b. It built the 1,848-panel domain and completed one frequency at `omega=1.2 rad/s`, producing current-source cache, mesh audit, and result files under the ignored example output area.

## 2. Before and after

| Item | Before | Release tree after cleanup |
|---|---:|---:|
| Git worktrees | 2 | 2 during this turn; the Phase 3b worktree is ready for removal |
| Local branches | 6 | 7 including the release branch; milestone branches remain explicit |
| Tracked files | 540-ish historical/production files | 125 |
| Non-ignored untracked files | 1,279 | 0 |
| Ignored/untracked runtime output in release tree | hundreds | generated example output only, ignored |
| Canonical examples | many historical cases | 1 single-sphere example with smoke/full configs |
| Absolute paths in active source/entry points | historical Phase 3 scripts contained them | 0 found in active source/runtime scan |

The first cleanup commit removed 415 tracked historical/generated files and added the portable entry points and release documentation. The research material was moved to `E:\CRESTU-1F_v1.0_20260824\Archive\CRESTU-ReleaseArchive-20260908` with a 2,039-file SHA-256 manifest.

## 3. Canonical baseline and Git history

The release is based on the complete `phase3-2b-response-convergence` chain, including Phase 3.2 and 3.2a ancestors, plus the audited checkpoint working changes. The key commits are:

```text
6e4740a  merge: preserve audited checkpoint working changes
f6f2c66  refactor: slim release tree and add portable entry points
ad219b3  fix: integrate portable first-order and reference parsing corrections
5cb2877  fix: add central numerical option defaults
```

The protected Phase 3b branch was pushed before cleanup:

```text
origin/phase3-2b-response-convergence -> 262892a
origin/codex/release-20260908        -> 5cb2877
```

The old checkpoint, Phase 3.2, Phase 3.2a branches and tags remain available for historical recovery. They are not required by a fresh clone.

## 4. What was integrated from CRESTU_GPT6PRO_FIXED

The source-aware review kept the independent, self-contained changes and rejected a blind package overlay.

Integrated production/source assets:

- `Source Code/1.Input/get_default_solver_options.m` — the missing central defaults helper required by the corrected panel-integral path.
- `Source Code/2.Mesh/audit_rankine_wave_resolution.m` — actual mesh/wavelength diagnostic.
- `Source Code/2.Mesh/validate_outer_waterline_compatibility.m` — explicit waterline seam guard.
- `Source Code/4.Potential/get_rankine_code_version.m` — conservative all-source cache fingerprinting.
- `Source Code/4.Potential/rankine_panel_integrals.m` — common planar treatment for warped quadrilateral single/double-layer integrals.
- `Source Code/5.Force/compute_haskind_excitation.m` — corrected diagnostic Haskind sign convention.
- `Source Code/5.Force/read_wamit_excitation.m`, `read_wamit_first_order.m`, `read_wamit_rao.m` — local DOF/ULEN/limit-record normalization.
- `Case_Wave/Common_Scripts/read_wamit_dataset.m` — portable raw-reference parsing and provenance handling.

Documentation integrated or rewritten:

- `AGENTS.md`
- `FIRST_ORDER_CONVENTIONS.md`
- `CRESTU_Theory_and_Technical_Manual.md` warning
- `docs/development/FIXED_PACKAGE_IMPORT_MANIFEST.csv`

Not integrated into production:

- The full `PARA14–16` parameterization sweep and its 30-file source change set. It is not yet a clean Phase 3b three-way merge and its delivery status still contains MATLAB-pending claims.
- Fixed-package full validation wrappers and large Python/MAT diagnostics. They remain in the local archive for later, explicit validation work.
- Large MAT/NPZ/CSV/PNG/WAMIT outputs and generated cache.
- Any content from `CRESTU_GPT6PRO_FIXED\GPT6PRO_ONE_SHOT\NextStep`. It was not read, copied, moved, uploaded, or modified.

## 5. Final release tree

```text
CRESTU/
├─ Source Code/                 production numerical source
├─ Case_Wave/                   five portable runtime/reference files
├─ Case_UniformFlow/            compact uniform-flow benchmark
├─ examples/single_sphere/      smoke/full config and body mesh
├─ validation/                  release validation guide
├─ docs/                        structure, conventions, known issues, import manifest
├─ scripts/maintenance/         maintenance policy
├─ AGENTS.md
├─ FIRST_ORDER_CONVENTIONS.md
├─ setup_crestu.m
├─ run_example_single_sphere.m
├─ run_validation_smoke.m
└─ README.md
```

Historical `Case_Wave` convergence, old cases, Phase 2/3 artifacts, shadows, handoff bundles, generated output, and cleanup reports are now outside the Git release tree in the local archive. Git history still preserves the committed Phase 3 chain and the pre-cleanup bundle preserves all refs.

## 6. Validation performed

### Repository/static checks

- Fresh clone from `origin/codex/release-20260908`: PASS.
- Fresh clone branch status: clean before generated example output.
- Fresh clone tracked file count: 125.
- `CRESTU_GPT6PRO_FIXED` and `NextStep` absent from the fresh clone.
- Active source/runtime absolute-path scan: no `E:\`, `C:\`, `CRESTU_GPT6PRO_FIXED`, or `NextStep` dependency found.

### MATLAB checks

- `setup_crestu`: PASS in the original release tree and fresh clone.
- `run_validation_smoke`: PASS; body 588 panels, free surface 280, bottom 868, outer 112, total 1,848 panels.
- `run_example_single_sphere('SMOKE')`: PASS in the original tree and fresh clone.
- Fresh-clone production assembly: PASS; 1,848×1,848 Rankine matrix, cache MISS, one frequency at `omega=1.2 rad/s`, result/cache/mesh-audit files created.
- Full 16-frequency Fine run: not executed in this release turn; it remains an explicit, more expensive command in README.

The one-frequency production run took approximately 279 s in the original tree and 343 s in the fresh clone. This is a runtime/caching/environment observation, not a performance benchmark.

## 7. Local archive

Archive location:

```text
E:\CRESTU-1F_v1.0_20260824\Archive\CRESTU-ReleaseArchive-20260908
```

The archive contains:

- moved full Case_Wave and old-case trees;
- Phase 2/3 artifacts and historical evidence;
- GPT6PRO_HANDOFF and the explicitly allowed GPT6PRO_ONE_SHOT package content;
- MATLAB preference/cache files;
- root cleanup/audit reports;
- `ARCHIVE_MANIFEST.csv` with relative path, size, timestamp and SHA-256.

The separate pre-cleanup bundle remains at:

```text
E:\CRESTU-1F_v1.0_20260824\Archive\PreCleanup-20260908\CRESTU_PRE_CLEANUP_20260908.bundle
```

## 8. How the junior should start

1. `git clone https://github.com/PengYunqiang/CRESTU.git`
2. Open MATLAB and `cd` to the cloned repository root.
3. Run `setup_crestu`.
4. Read `AGENTS.md` and `FIRST_ORDER_CONVENTIONS.md`.
5. Run `report = run_validation_smoke`.
6. Run `results = run_example_single_sphere` for the one-frequency production smoke.
7. Inspect generated files under `examples/single_sphere/`; they are ignored and local.
8. Run `results = run_example_single_sphere('FULL')` when the smoke case is understood.
9. Put new small reproducible inputs under `examples/`, runners under `validation/`, and generated outputs outside tracked source.
10. Record numerical changes, source/config identity, and validation status before committing.

## 9. Remaining release boundary

The release is handoff-ready for source understanding, setup, modification, and a representative runtime. It does not claim that every historical WAMIT comparison, Phase 3 convergence result, or external fixed-package MATLAB validation is complete. Those are deliberately separated into the local archive and documented as partial/pending rather than hidden in the release.
