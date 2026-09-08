# Release validation

The release validation surface is intentionally small:

- `run_validation_smoke.m` checks setup and canonical geometry.
- `run_example_single_sphere.m` runs the production solver on one frequency.
- `examples/single_sphere/Case1_Fine.cfg` provides the full canonical 16-frequency input.

Large Phase 2/3 runners and generated evidence are retained in Git history or the local archive, not as fresh-clone dependencies. Add a new validation runner only when it has a reproducible input, a declared pass criterion, and no dependency on historical cache.
