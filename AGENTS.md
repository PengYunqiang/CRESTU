# CRESTU Project Rules

## Purpose

CRESTU is MATLAB research software implementing a three-dimensional Rankine-source BEM hydrodynamic solver. The repository must remain cloneable, reproducible, and understandable without access to the original developer machine.

## Canonical runtime

- Start from the repository root and run `setup_crestu` before calling project functions.
- The production entry point is `Case_Wave/run_frequency_domain_case.m`.
- The beginner entry point is `run_example_single_sphere`.
- Canonical source remains under `Source Code/`; do not add personal absolute paths to MATLAB files or configuration files.
- Generated results belong under ignored example/output directories, not beside production source.

## Physical conventions

- SI units throughout.
- Harmonic convention: `exp(+i*omega*t)`.
- Body normals point outward from the body and into the fluid: `n_B`.
- DOF order per body: Surge, Sway, Heave, Roll, Pitch, Yaw.
- Generalized normals are `N=[n_B,(x-CG) cross n_B]`.
- Radiation potential is per unit generalized displacement; radiation RHS is `-S*(i*omega*N)`.
- Production radiation damping is the pressure-integrated path. Haskind and control-cylinder flux are diagnostics, not replacements.
- Read `FIRST_ORDER_CONVENTIONS.md` before changing a sign, normal, phase, or A/B extraction formula.

## Data and validation

- `examples/` contains small, canonical, reproducible inputs only.
- `validation/` contains runners, fixtures, and validation documentation; large results remain external.
- WAMIT raw data is optional external reference data unless a small fixture is explicitly committed under `validation/fixtures/`.
- Historical reports and MAT/CSV results are forensic evidence, not proof of current source behavior.
- Never fabricate a missing body mesh or silently fall back to a different geometry.
- Do not commit PotCache, RunRecord, large MAT, solver output, logs, or generated meshes unless a fixture is explicitly required for a test.
- Preserve complex force/RAO phase, source/config provenance, geometry counts, and code fingerprints when adding validation.

## Editing and commits

- Keep production changes minimal and physically motivated.
- Separate source fixes, validation changes, documentation, and repository cleanup into logical commits.
- Do not use `git add .` for research cleanup; stage explicit paths.
- Do not modify `CRESTU_GPT6PRO_FIXED/GPT6PRO_ONE_SHOT/NextStep` or make the release depend on it. That directory is future work and is outside this repository's release boundary.
- If a numerical result is not verified under the current MATLAB source, label it `PENDING`, `PARTIAL`, or `UNKNOWN`; do not call it validated.

## Release expectation

A fresh clone must be able to run `setup_crestu`, execute `run_example_single_sphere`, and run `run_validation_smoke` without developer-specific absolute paths. If MATLAB runtime validation is unavailable, keep the static/path checks and document the pending runtime explicitly.
