# Known issues and validation boundaries

This file separates release packaging decisions from unresolved scientific questions.

## MATLAB runtime status

The fixed GPT6PRO delivery contains strong independent Python evidence and prepared MATLAB runners, but its own delivery checks state that MATLAB execution was not performed. Do not label those Python metrics as MATLAB validation. A real MATLAB smoke/full run is the release gate for numerical claims.

## WAMIT reference geometry

The supplied FullSphereIRR3 reference uses a 300-panel hemisphere. The canonical Fine body input uses 588 panels. This is a documented geometry difference, not an automatic solver defect. Reference comparison must record the exact body mesh, input deck, frequency conversion, and source hashes.

## Phase 3 mesh/response studies

The Phase 3.2 → 3.2a → 3.2b chain is preserved in Git and remote branches. Its formal convergence and response gates remain research evidence, not a blanket claim that every geometry or frequency is physically validated.

## Outer-domain BMF history

Historical OuterDomain Fine BMF files had substantial frequency-local geometry changes. They are not release inputs and should not be copied over canonical example geometry without provenance.

## Generated results

Large MAT, PotCache, RunRecord, solver output, logs, and sweep tables are intentionally excluded from the release. Recreate them with a documented configuration and record the resulting source/config/reference identity.
