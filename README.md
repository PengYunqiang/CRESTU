# CRESTU-1F

CRESTU-1F is a MATLAB Rankine-source boundary-element solver for research hydrodynamics. It supports single- and multi-body frequency-domain radiation/diffraction calculations, added mass, pressure-integrated radiation damping, wave excitation, hydrostatics, RAOs, symmetry reduction, and selected validation workflows.

The repository is an active research release. Numerical results must be interpreted together with the validation status and limitations documented under `docs/` and `validation/`.

## Quick start

Requirements:

- MATLAB R2022b or newer; the maintained development environment uses MATLAB R2023b.
- No external data is required for the minimal smoke example.
- WAMIT is optional and only needed for reference comparisons; see `validation/VALIDATION_GUIDE.md`.

From MATLAB, starting in the repository root:

```matlab
setup_crestu
results = run_example_single_sphere
```

The example uses the canonical single-sphere geometry in `examples/single_sphere/`. Generated meshes, cache, results, and logs in that directory are ignored by Git.

For a geometry/path validation smoke test:

```matlab
run_validation_smoke
```

For the full 16-frequency Fine configuration, use the explicit configuration after the smoke run has succeeded:

```matlab
setup_crestu
cfg = fullfile('examples','single_sphere','Case1_Fine.cfg');
results = run_frequency_domain_case(cfg, skipPhysicalDiagnostics=true);
```

## Repository structure

```text
Source Code/                 production MATLAB source, grouped by numerical role
Case_Wave/                   historical and maintained wave-case runners
Case_UniformFlow/            uniform-flow benchmark code
examples/single_sphere/      small canonical input and smoke case
validation/                  validation guidance, runners, fixtures, provenance
docs/                        structure, conventions, known issues, development notes
scripts/                     maintenance and release helpers
FIRST_ORDER_CONVENTIONS.md   normative first-order equations and signs
setup_crestu.m               portable path setup
run_example_single_sphere.m  beginner-facing example entry point
run_validation_smoke.m       fast geometry/configuration smoke test
```

## Physical conventions

The normative conventions are in [FIRST_ORDER_CONVENTIONS.md](FIRST_ORDER_CONVENTIONS.md). In short: SI units, `exp(+i*omega*t)`, body-outward normals, DOFs ordered Surge/Sway/Heave/Roll/Pitch/Yaw, and unit-displacement radiation. Production damping is pressure-integrated; Haskind and control-cylinder flux paths remain diagnostics.

## Validation and external reference data

Start with [validation/VALIDATION_GUIDE.md](validation/VALIDATION_GUIDE.md). The repository contains small canonical inputs and runners, not the entire research output history. Full WAMIT data is an optional external dataset. Expected folder layout, file roles, and provenance requirements are documented there.

The external `CRESTU_GPT6PRO_FIXED` delivery, historical Phase 2/3 evidence, and large generated outputs are local archive material. They are not runtime dependencies of a fresh clone. In particular, `NextStep` is future development and is intentionally outside the release boundary.

## Development

1. Read [AGENTS.md](AGENTS.md) and [docs/REPOSITORY_STRUCTURE.md](docs/REPOSITORY_STRUCTURE.md).
2. Run `setup_crestu` from the repository root.
3. Add new small, reproducible inputs under `examples/` and validation runners under `validation/runners/`.
4. Keep generated output under ignored output/cache directories.
5. Preserve source/config/reference provenance and complex force/RAO phase in numerical validation.
6. Record unresolved scientific questions in [docs/development/KNOWN_ISSUES.md](docs/development/KNOWN_ISSUES.md).

## License

CRESTU is distributed under GPLv3. See [LICENSE](LICENSE).
