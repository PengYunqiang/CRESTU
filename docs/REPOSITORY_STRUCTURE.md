# Repository structure

CRESTU intentionally keeps its MATLAB production paths stable. The release is organized around discoverability and portability without renaming the existing numerical source directories.

```text
CRESTU/
├─ Source Code/                 production numerical source
│  ├─ 1.Input/                  config parsing and validation
│  ├─ 2.Mesh/                   BMF geometry and boundary generation
│  ├─ 3.HessSmith/              steady source-panel kernels
│  ├─ 4.Potential/              Rankine BIE, dispersion, RHS, diagnostics
│  ├─ 5.Force/                  A/B, excitation, hydrostatics, RAO
│  └─ 6.MeanDriftLoads/         second-order drift utilities
├─ Case_Wave/                   maintained wave runners and historical cases
├─ Case_UniformFlow/            uniform-flow benchmark
├─ examples/single_sphere/      canonical small input and smoke case
├─ validation/                  validation guide and release-facing runners
├─ docs/                        conventions, structure, known issues, history
├─ scripts/                     maintenance helpers
├─ setup_crestu.m               portable path setup
├─ run_example_single_sphere.m  beginner-facing runtime entry point
└─ run_validation_smoke.m       fast geometry/configuration check
```

## Runtime paths

`setup_crestu` resolves the repository root from its own file location and adds the required source directories. Production code should use configuration-relative paths and `fullfile`; it must not depend on the caller's current directory or a developer-specific absolute path.

## Data policy

- Tracked: source, small canonical inputs, reproducible runners, conventions, validation metadata, and concise reports.
- Ignored: generated BMF, cache, MAT results, logs, and output directories.
- External: full WAMIT archives, historical Phase 2/3 evidence, large research results, and GPT6PRO delivery packages.
- Excluded from this release: `CRESTU_GPT6PRO_FIXED/GPT6PRO_ONE_SHOT/NextStep`.
