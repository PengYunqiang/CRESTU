# Validation guide

## Fast checks

From the repository root in MATLAB:

```matlab
setup_crestu
report = run_validation_smoke
results = run_example_single_sphere
```

`run_validation_smoke` builds the canonical one-frequency domain and checks that body, free-surface, bottom, and outer-boundary panels exist. `run_example_single_sphere` then runs the one-frequency production path with physical diagnostics skipped so that a clone can establish a basic runtime quickly.

## Full canonical example

After the smoke run succeeds:

```matlab
results = run_example_single_sphere('FULL')
```

This uses the 16-frequency `examples/single_sphere/Case1_Fine.cfg`. It is more expensive than the smoke case and writes generated output beside the example configuration; those outputs are ignored by Git.

## WAMIT comparison

WAMIT is optional for the basic example. For reference comparisons, provide a local external dataset with this structure:

```text
WAMIT_REFERENCE/
├─ FullSphereIRR3/
│  ├─ hemisphere.1
│  ├─ hemisphere.2
│  ├─ hemisphere.4
│  ├─ hemiSphere.gdf
│  ├─ hemiSphere.POT
│  └─ hemiSphere.frc
└─ README.md
```

The full local WAMIT archive is outside Git. Record its path, role, geometry, units, period-to-frequency conversion, and file identity in a provenance document. Do not silently substitute a different reference.

## Evidence levels

- `RUNTIME_TESTED`: a real MATLAB run completed.
- `NUMERICALLY_VERIFIED`: a representative result passed a declared regression/analytic metric.
- `PHYSICALLY_VALIDATED`: reserved for evidence against a trusted benchmark or experiment with matching geometry and conventions.
- `PENDING`/`PARTIAL`: prepared evidence or auxiliary Python evidence that still needs the current MATLAB path.

## Historical validation

Phase 2/3 evidence remains available in the local archive and preserved Git history. It is not copied into the minimal example or required for a fresh clone.
