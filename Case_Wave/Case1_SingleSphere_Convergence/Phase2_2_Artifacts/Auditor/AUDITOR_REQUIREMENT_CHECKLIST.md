# Phase 2.2 Auditor Requirement Checklist

> **Artifact class:** independent audit metadata; not production code, solver configuration, test implementation, or Executor output.
>
> **Owner:** Agent B — Auditor. Production-code edits are prohibited for this owner.
>
> **Status legend:** `PENDING`, `PASS`, `WARNING`, `BLOCKER`, `N/A`. A checkpoint cannot pass while any in-scope item is `BLOCKER`. Every blocker must be cleared by Executor evidence before final acceptance.

## 0. Audit control

| ID | Requirement / acceptance evidence | Status | Evidence / note |
|---|---|---|---|
| AC-01 | Preserve the confirmed Phase 2.1 fixes: separated body/outer studies; geometry-bound result/potential cache keys; raw/projected/energy separation; RAO not overwritten by energy/Haskind; clean/reload machine agreement. Do not reopen or alter them without new evidence. | PASS | Previous report read; Executor explicitly freezes these areas. Any later diff touching them reopens this item. |
| AC-02 | Scope is single-sphere, first-order Rankine Phase 2.2. No multi-body, roll-specialist, second-order drift/QTF work. | PASS | Executor plan is confined to the requested single-sphere first-order scope. |
| AC-03 | Executor alone may edit production code. Auditor may update only independent audit metadata/review. | PASS | This file is isolated under `Phase2_2_Artifacts/Auditor/`. |
| AC-04 | Auditor reviews requirements, plan, critical diffs/formulae/configurations, decisive numerical evidence, test validity, final report; it does not duplicate full-repo review or large computations. | PASS | Review protocol established below. |
| AC-05 | Conclusions must identify unresolved blockers honestly; `PARTIAL` is required where evidence is partial. Phase may be declared complete only with `FINAL_AUDIT = PASS`. | PENDING | |

## 1. Baseline and experimental controls

| ID | Requirement / acceptance evidence | Status | Evidence / note |
|---|---|---|---|
| BL-01 | Read current project/tests/latest `NUMERICAL_AUDIT_REPORT.md`; continue from confirmed results rather than restarting or refactoring validated areas. | PASS | Auditor read latest report. Executor reports reading the report plus relevant single-sphere solver/mesh/BC/force/cache/test call chain and continues from confirmed fixes. |
| BL-02 | Fix the current Fine body mesh throughout Phase 2.2 diagnostics. Do not change body geometry, body-pressure integration, basic RAO definition, verified cache fingerprint, or raw/projected separation. | PENDING | Plan fixes Fine body at 588 panels and freezes protected logic; verify actual body count/hash and absence of protected diffs in Checkpoints 2/4. |
| BL-03 | Baseline includes at least one low-, one mid-, and at least three representative `omega >= 1.5 rad/s` frequencies, with actual values stated. | PENDING | Plan specifies omega = {0.5, 1.1, 1.5, 1.7, 2.0}; actual clean-run evidence pending. |
| BL-04 | Every decisive diagnostic uses clean recompute and records cache miss; old cache is not used to infer root cause. Cache reload may only be a separate reproducibility check. | PENDING | Plan specifies isolated output and `IPOTEN=1`; Checkpoint 2 must prove this setting causes actual clean recompute/cache MISS rather than merely requesting potential output. |
| BL-05 | Save a pre-change baseline before each evidence-supported production fix and retain before/after comparable data. | PENDING | Explicitly included in plan; artifacts pending. |

## 2. Outer-domain OFAT and geometry authenticity

| ID | Requirement / acceptance evidence | Status | Evidence / note |
|---|---|---|---|
| OD-01 | A genuine OFAT design varies exactly one of: FS panel spacing; FS radial extent; bottom panel spacing; bottom radial extent; outer-vertical-boundary spacing; outer radius; sponge/radiation-zone start; zone width; current sponge/radiation parameters. Other factors are demonstrably fixed. | PENDING | Revised design covers all nine factors. FS/bottom physical radii use independent top/bottom frustum controls; connector geometry/tessellation is a deterministic closure consequence and results must be named top/bottom truncation-shape sensitivity. Frustum BC must pass Checkpoint 3 before those two studies run or count. Actual invariant/delta manifests remain pending. |
| OD-02 | Each run records body/FS/bottom/outer panel counts, total unknowns, component/whole mesh SHA/hash, cache hit/miss, frequency, wavelength, and matrix dimensions. | PENDING | All fields are in the plan; actual rows pending. Whole hash must include BC-affecting non-geometric parameters or be accompanied by a configuration/BC hash. |
| OD-03 | Configuration changes are proven to alter the intended active solver geometry/BC only; no fake experiment from unchanged mesh, wrong object, stale cache, rounded-away setting, or wrong execution path. | PENDING | Design includes isolated new paths, nonexistence/preexisting-cache assertions, observed MISS, component/merged/config/code/BC/cache SHA, canonical effective manifests, allowed-field deltas, and expected hash changes. Actual manifests/rows remain required. For factors 7–9, geometry hashes must remain identical while effective BC/damping and cache hashes change. |
| OD-04 | Radius and resolution are not conflated into a generic Coarse/Medium/Fine convergence claim. | PENDING | Orthogonal q/alpha design separates common physical truncation radius from component resolution. However the two proposed “resolved radial extent” studies are grading studies and must remain clearly classified under resolution, not radius. |
| OD-05 | Radius × Resolution A: hold `h/lambda` approximately constant while increasing outer radius/lambda; quantify how closely resolution is held. | PENDING | Design: alpha=0.25, q={1,1.5,2,2.5}; actual componentwise resolution tolerance <=5%. Executed effective geometry/spacing evidence pending. |
| OD-06 | Radius × Resolution B: hold geometric extent/radius fixed while changing `h/lambda`; prove identical radii/hashes except intended remeshing. | PENDING | Design: q=1.5, alpha={0.35,0.25,0.18,0.125}. Same physical component extents must be shown from manifests; remeshing hashes are expected to differ. |
| OD-07 | Radius × Resolution C: improve radius and resolution together by an explicit stated ratio; interpret only as mixed-error evidence. | PENDING | Design pairs={(1,0.35),(1.5,0.25),(2,0.18),(2.5,0.125)}; correctly classified as mixed-error evidence. Execution pending. |
| OD-08 | Do not call movement toward WAMIT “convergence”; convergence must be supported by successive self-consistent numerical refinements and asymptotic behavior/error evidence. | PENDING | |

## 3. Wavelength-normalized resolution

| ID | Requirement / acceptance evidence | Status | Evidence / note |
|---|---|---|---|
| WR-01 | Dispersion relation and computed wavelength are stated (including finite depth as applicable), traceable, and dimensionally consistent. | PENDING | Finite-depth dispersion explicitly planned; derivation/source and results pending. |
| WR-02 | Report `h/L_body`, `h/lambda`, `lambda/h` (panels per wavelength), `R_outer/lambda`, `R_FS/lambda`, and `R_bottom/lambda` for relevant components/runs. | PENDING | Design reports component/band radial, tangential and vertical scales plus P95/max/aspect/skew rather than an aggregate only. Results pending. |
| WR-03 | Explicitly test whether nominal “Fine” becomes under-resolved as omega rises and lambda falls. | PENDING | Design audit estimates worst directional high-frequency resolution at only about 1.23–3.18 panels/lambda, a credible under-resolution hypothesis. It is not yet a causal conclusion. |
| WR-04 | Deliver high-frequency `quantity vs h/lambda` and `error/sensitivity vs h/lambda` data and figures, using unsmoothed raw physical values. | PENDING | |

## 4. Boundary-condition mathematical audit

| ID | Requirement / acceptance evidence | Status | Evidence / note |
|---|---|---|---|
| BC-01 | Declare whether code uses `exp(+i omega t)` or `exp(-i omega t)` and use one convention consistently in all derivations and labels. | PASS | Incident `exp(-ikx)` propagates in +x under `exp(+iωt)`; `Phi_I=i g a/ω V(z) exp(-ikx)` and `p=-iρ ω Phi` are consistent for peak amplitude `a`. |
| BC-02 | For radiation body Neumann, diffraction body Neumann, FS Robin, bottom Neumann, outer outgoing-wave BC, radiation potential→pressure→A/B, and diffraction potential→excitation: provide theory → time convention → normal definition → exact code location mapping. | BLOCKER | Mapping exposes production defects: radiation RHS sign; FS source orientation/finite-depth coefficient; outer BIE sign/radial term/sponge-consistent wavenumber; pressure-to-force sign. Diffraction body RHS and bottom zero-flux mapping are consistent. Focused tests and fixes are required before numerical studies. |
| BC-03 | Independently trace `+iω/-iω`, incident phase, unit-displacement versus unit-velocity radiation normalization, wave-amplitude normalization, diffraction RHS sign, and real/imag mapping to A/B. | BLOCKER | Time/incident/peak normalization and diffraction RHS pass analytically. Radiation columns are unit generalized displacement (`q_s=iω n_j`); correct all-inward-source BIE requires `rhs=-S iω n_j`. Physical body force is `Q=+iρ ω integral(Phi n_j dS)`. Current opposite RHS and opposite force signs can cancel in A/B and must be fixed/tested together. |
| BC-04 | Define body, source-panel, collocation, and outer-boundary normal orientations and verify implementation signs/self terms against those definitions. | BLOCKER | `n_B` points body→fluid and equals `-n_Omega`; bottom +ez and outer -er are also `-n_Omega`, but FS +ez equals `+n_Omega`. Code `D_s=1/2-K_s` assumes every source normal is `-n_Omega`. FS is inconsistent. Kernel derivative also depends on vertex winding, so any normal flip must flip/validate winding, not metadata alone. Collocation has no normal. |
| BC-05 | Audit Sommerfeld `ik` sign, radial-decay term, finite-depth dispersion/depth factor, and sponge/radiation-zone implementation. | BLOCKER | Physical FS coefficient is `nu0=ω²/g=k tanh(kh)`, not real `k` at finite depth. For cylinder with inward source normal and `exp(+iωt)`, outgoing `exp(-iκr)/sqrt(r)` gives `q_s/Phi=iκ+1/(2r)` and BIE row `D_s+S q_s=0`; current row has opposite sign and omits decay. Frustum adds `R' κ tanh[κ(z+h)]` divided by `sqrt(1+R'^2)`. In sponge plateau, κ must be the complex root consistent with the chosen complex FS Robin coefficient; using the undamped real physical k is not a closed derivation. |
| BC-06 | Any convention-mixing claim is first demonstrated by a minimal reproduction/unit test before production modification. | BLOCKER | T1–T11 are proposed but not run. The prioritized mandatory sequence is recorded under Checkpoint 3. No production modification is authorized yet. |

## 5. Three independent damping closures and dimensional audit

| ID | Requirement / acceptance evidence | Status | Evidence / note |
|---|---|---|---|
| DP-01 | Persist independent per-frequency `B_pressure` from radiation body-pressure integration; it is not overwritten by other paths. | PENDING | Planned; existing raw pressure path must be shown not to alias later closures. |
| DP-02 | Persist independent per-frequency `B_flux` from radiated-wave/far-field energy flux; it does not alias pressure or Haskind data. | PENDING | Planned as a direct far-field flux. Executor correctly flags that the existing “energy” value is actually a Haskind outer product, so it cannot be accepted as `B_flux`. |
| DP-03 | Persist independent per-frequency `B_haskind` only when theory and normalization are confirmed; otherwise explicitly label `NOT VALIDATED` (or final `NOT_APPLICABLE`) without substitution. | PENDING | Plan preserves an independent Haskind path or reports `NOT VALIDATED`; implementation and normalization proof pending. |
| DP-04 | Output `B_pressure`, `B_flux`, `B_haskind`, and all three pairwise ratios by frequency, with invalid/zero denominators handled transparently. | PENDING | Flux must retain the raw cross-mode pressure-work bilinear before any real/symmetric presentation; no clipping/projection may hide it. |
| DP-05 | Audit each damping prefactor/normalization term: `rho`, `g`, `omega`, `k`, phase/group velocity, finite-depth factor, `tanh(kh)`, amplitude basis, displacement/velocity basis, peak/RMS, complex convention, factors 2/pi/2pi, angular integration. | BLOCKER | Diagonal control-surface identity is `P_out=1/2 Re integral(p q_out*)`, `p=-iρ ω Phi`, and for unit displacement `P_out=1/2 ω² B`, hence `B=(ρ/ω) Im integral(Phi q_out*)`. Cross-mode B requires the correctly indexed Hermitian/real bilinear, not blind elementwise `imag`. Group-velocity plane-wave test must fix the peak-amplitude factor 1/2. Haskind prefactor remains unvalidated. |
| DP-06 | Provide a B33 dimensional derivation/check; numerical agreement alone is not validation. | PENDING | Required derivation must include division by generalized displacement amplitudes: `A=Re(Q)/(ω² ξ)`, `B=-Im(Q)/(ω ξ)` and `B_flux=Re(C)/(ω² ξ_i ξ_j*)`; solver currently assumes numerical 1 m/1 rad without unit-carrying interfaces. |
| DP-07 | Determine which closure first loses agreement and separate a faulty diagnostic formula from a faulty solver BC/solution. | BLOCKER | Existing `compute_radiation_damping_energy` is a Haskind-force outer product, not direct flux. Current Haskind body integrand uses a plus where the same-normal Green concomitant is antisymmetric (`Phi_I q_R - Phi_R q_I`); it is `NOT VALIDATED` until a manufactured identity fixes sign, conjugation, incident direction and prefactor. |

## 6. Raw physics diagnostics and conditioning

| ID | Requirement / acceptance evidence | Status | Evidence / note |
|---|---|---|---|
| RD-01 | For every representative frequency save `A_raw`, `B_raw`, raw reciprocity residual, diagonal B, symmetric B, and `min eig(sym(B_raw))`; projected values do not hide raw behavior. | PENDING | |
| RD-02 | Save radiation and diffraction linear residuals, solution norm, maximum radiation potential, maximum FS potential, and outer-boundary potential magnitude. | PENDING | |
| RD-03 | Save BC residual statistics separately for body, FS, bottom, and outer boundary; define norm/statistics and normalization. No global-only residual. | PENDING | Component split is planned. Checkpoint 3 must distinguish algebraic row residual from independently evaluated physical BC residual; a solved-row residual alone is circular and insufficient for local physics validation. |
| RD-04 | Record a reliable matrix condition estimate if algorithmically feasible; if infeasible, state method limitation and use clearly labelled proxy evidence rather than inventing a condition number. | PENDING | Reliable estimate is planned; method, scaling sensitivity, estimator meaning, and feasibility pending. |
| RD-05 | Correlate conditioning/residual/solution/potential diagnostics with panel count, radius, `h/lambda`, and B33/F3/RAO3 sensitivity. | PENDING | |
| RD-06 | If refinement worsens results, discriminate among non-asymptotic mesh, conditioning, BC error, quadrature/panel-integral error, wave under-resolution, and boundary reflection; do not assume Medium or Fine is correct. | PENDING | |

## 7. WAMIT validation and forbidden shortcuts

| ID | Requirement / acceptance evidence | Status | Evidence / note |
|---|---|---|---|
| WV-01 | WAMIT is independent validation only, never a tuning/selection target. CRESTU and WAMIT remain unscaled raw physical values. | PENDING | |
| WV-02 | Each comparison includes A33, B33, `|F3|`, F3 phase, RAO amplitude, and RAO phase where available—not amplitude alone. | PENDING | |
| FS-01 | No frequency-dependent empirical correction, WAMIT-based scale, artificial multiplier, smoothing, curve fit, damping tuned for RAO, or plotting-data manipulation. | PENDING | |
| FS-02 | No reciprocity projection to conceal raw-solver defects; no PSD clipping to conceal negative damping. | PENDING | |
| FS-03 | No prohibited scope diversion or body-mesh refinement as the proposed cure. | PENDING | |

## 8. Production-change evidence gate

| ID | Requirement / acceptance evidence | Status | Evidence / note |
|---|---|---|---|
| CH-01 | Forensic audit precedes any production change; every change identifies one specific physical/numerical defect with theoretical or direct numerical evidence. | PASS | First proposed change is supported by Green derivation and fail-before T2: all-inward cube residual `5.55e-17`, current mixed-orientation residual `0.21795`, with winding alignment controls passing. Authorization is limited to double-layer source-orientation canonicalization. |
| CH-02 | One physical issue is changed at a time; no bundled unrelated fixes that prevent attribution. | PENDING | Proposed first patch changes only FS double-layer source orientation plus mandatory metadata/cache invalidation. FS Robin, radiation RHS, force, outer BC and damping closures must remain unchanged in this patch. Verify diff after implementation. |
| CH-03 | Each change records exact code location, why old implementation is wrong, governing formula, and a focused valid unit/minimal test that would fail before and pass after. | PENDING | Fail-before evidence is valid and preserved. After-test must exercise production canonicalization (shared routine or production-exposed orientation metadata plus inspected assembly path), not merely duplicate the multiplication in a test helper. |
| CH-04 | Each change reruns representative clean frequencies and only the necessary convergence studies; presents before/after evidence and receives Auditor review. | PENDING | |
| CH-05 | Test validity is audited: correct target/call path/config, assertions sensitive to defect, no cache contamination, intended geometry actually changed, and no circular oracle. | PENDING | |

## 9. Deliverables and decision questions

| ID | Requirement / acceptance evidence | Status | Evidence / note |
|---|---|---|---|
| DL-01 | Deliver `PHASE2_2_OUTER_DOMAIN_AUDIT.md` with all 18 requested sections, including this requirement checklist and an independent Auditor Review. | PENDING | |
| DL-02 | Retain supporting CSV, MAT, convergence figures, BC residual diagnostics, damping-closure comparison, and before/after comparison; report links/paths resolve. | PENDING | |
| DL-03 | Report directly answers Q1 root cause of ~96% high-frequency B33 outer sensitivity and Q2 apportions radius/resolution/radiation BC/FS BC/extraction/flux/conditioning or mixed causes with data. | PENDING | |
| DL-04 | Report answers Q3 first failing damping path, Q4 whether F3 shares B33 root cause, Q5 whether Fine outer is asymptotically converged, and Q6 whether ready for multi-body. | PENDING | |
| DL-05 | Executor prints exactly the nine requested final state variables with `YES/NO/PARTIAL/NOT_APPLICABLE` only as allowed by each field. | PENDING | |
| DL-06 | Auditor independently prints: `REQUIREMENTS_COMPLETE`, `PHYSICS_CHECK`, `NUMERICAL_VALIDATION`, `CONVERGENCE_CLAIM_JUSTIFIED`, `FORBIDDEN_SHORTCUTS_USED`, `UNRESOLVED_BLOCKERS`, `FINAL_AUDIT`. | PENDING | |
| DL-07 | Final baseline is reproducible, explainable, and demonstrates genuine numerical convergence where claimed; no premature multi-body transition. | PENDING | |

## 10. Mandatory checkpoint protocol

### Checkpoint 1 — Executor plan

Evidence requested: a concise plan naming baseline frequencies, clean-recompute controls, nine-factor OFAT matrix, radius×resolution design, BC derivation/damping-closure workflow, raw diagnostics, expected artifacts, change gate, and review handoffs.

Decision rule: `PASS` only if the plan covers the checklist without prohibited scope or premature production edits. Gaps that can be corrected before execution are `BLOCKER`, not deferred silently.

**Current verdict:** `PASS`.

Plan coverage is complete enough to proceed. No production change or numerical conclusion is authorized by this verdict. Execution evidence remains pending in the tables above.

Carry-forward warnings (must be resolved at the named later checkpoints):

1. `IPOTEN=1` is not by itself proof of a clean recompute; Checkpoint 2 requires observed cache MISS, isolated outputs, and keys/hashes tied to every varied geometry/BC parameter.
2. The old `PARA9` coupling, frequency-dependent automatic outer expansion, and dynamic `mu` can silently violate OFAT. Checkpoint 2 requires an explicit invariant/delta table of the *effective* per-frequency parameters, not only input configs.
3. Exact radius/resolution levels and the tolerance used to call `h/lambda` constant are deferred to Checkpoint 2.
4. Component BC residuals must not merely partition the solved matrix residual. Checkpoint 3 requires a definition showing what independent physical defect they can detect; otherwise they must be labelled algebraic residuals only.
5. Existing “energy” output is not an independent flux closure. `B_flux` must be newly derived/evaluated from actual energy flux, while Haskind remains separate or `NOT VALIDATED`.
6. A condition estimate must state matrix scaling and estimator limitations; if it is not reliable, report unavailability/proxies rather than a spurious precise number.

### Checkpoint 2 — Outer-domain design

Audit OFAT isolation, actual active-mesh change, full per-run provenance, wavelength normalization, three orthogonal groups, and clean-cache discipline.

**Current verdict:** `PASS`.

What passes at design level:

- Fine body, five frequencies, frequency-isolated clean paths, explicit effective manifests and SHA/cache assertions are well controlled. `IPOTEN=1` is supplemented by hard absence/preexisting-cache checks and observed `cacheStatus=miss`.
- Automatic expansion and dynamic `mu` are explicitly frozen. Non-geometric factors carry BC/effective/cache hashes even when mesh hashes must remain identical.
- FS spacing is controlled by one scalar `s_FS` with a deterministic 56-node transition, and will report directional/bandwise `h/lambda` and quality statistics. Closure and transition-mesh quality remain execution gates.
- Sponge start and width are now independently parameterized: a zero region, quadratic ramp of fixed start or width, and a plateau to fixed `R_outer`. This resolves the former identity `width=R_outer-r_inner`. Factor 7 fixes width/Router/mu; factor 8 fixes start/Router/mu; geometry must be bit-identical.
- Radius×Resolution A/B/C have explicit q/alpha levels, include FS radial/tangential, bottom, and outer vertical/tangential resolution, and impose <=5% componentwise tolerance rather than an aggregate scale.
- Algebraic component residuals are honestly labelled; the proposed off-surface/Richardson residual is separate and may be marked `NOT VALIDATED` if it fails. Flux and Haskind paths remain separated. Raw/scaled `rcond` will carry limitations.
- FS and bottom physical radial extents are now independent controls on a watertight frustum: `R(z)=R_b+(z+h)(R_t-R_b)/h`. Factor 2 changes only `R_t=20+q_t lambda` with `q_t={1.5,2,2.5}` and fixed `R_b=20+1.5 lambda`; factor 4 is symmetric. Wall shape/normal/area/count/hash are explicitly derived closure effects, not extra independent inputs. Conclusions are restricted to top/bottom truncation-shape sensitivity.
- Connector resolution is no longer a hidden confound. Its meridional count is derived from `L_s=sqrt(h^2+(R_t-R_b)^2)` and a fixed 12.5 m target; tangential count is derived from a fixed baseline arc target and rounded only toward finer multiples of four. Actual top/bottom arc/chord and meridional `h/lambda`, exceptions, and <=5% acceptance are required. Factor 5 separately changes wall spacing on a fixed cylinder.
- Watertight shared edges, scale-normalized ring mismatch, Euler/connectivity, frustum area/volume, analytic source-normal, winding/orientation, and self-intersection assertions are mandatory before assembly.

Carry-forward gate: the frustum FARFIELD outgoing condition must be derived for its nonzero vertical normal component under the verified time and source-normal conventions, then pass Checkpoint 3 and a focused test. Until that happens, factors 2 and 4 must not run and cannot support attribution. A failure here reopens Checkpoint 2 as a blocker rather than permitting the cylindrical-wall formula to be reused.

### Checkpoint 3 — BC mathematics

Independently check time convention, `iω`, normals, radiation/diffraction/FS/bottom/outer BCs, pressure and A/B extraction, excitation, amplitude and finite-depth normalization.

**Current verdict:** `BLOCKER`.

Independent derivation and findings:

1. With `G=1/r` and fluid-domain outward normal `n_Omega`, Green's boundary equation is
   `(1/2 I + K_Omega) Phi = S q_Omega`.
   If every stored/source normal is `n_s=-n_Omega`, then `K_s=-K_Omega`, `q_s=-q_Omega`, and the equivalent equation is
   `(1/2 I-K_s) Phi = -S q_s`.
   The code's `D_s=1/2-K_s` therefore does **not** permit arbitrary mixed source orientations.
2. The body, bottom and cylindrical outer source normals are inward to the fluid domain (`-n_Omega`); FS `+e_z` is outward from the fluid (`+n_Omega`). Thus only the FS violates the matrix's assumed source orientation. Because `rankine_panel_integrals` obtains the solid-angle sign from vertex order, changing only the stored normal is insufficient; supplied normal, winding and analytic derivative must agree.
3. Radiation for unit generalized displacement has body/source derivative `q_R,s=iω n_j`. The correct RHS is `-S iω n_j`; the current `+S iω n_j` is reversed.
4. Fixed-body diffraction has `q_D,s=-q_I,s`; substitution gives `rhs=+S q_I,s`, so the current diffraction body RHS is correct under the stated body normal.
5. On the physical FS, `partial_z Phi=(ω²/g)Phi=k tanh(kh)Phi`. The real part of the current `nu=k` is not the finite-depth coefficient except asymptotically. A sponge may make the Robin coefficient complex, but its exact PDE and sign must be stated; it cannot remove the physical-region factor `tanh(kh)`.
6. Bottom `q_s=0` with stored `+e_z=-n_Omega` is consistent once the global source convention is enforced.
7. For `exp(+iωt)`, a cylindrical outgoing mode is `Phi~f_kappa(z) exp(-i kappa r)/sqrt(r)`. With inward source normal, `q_s/Phi=i kappa+1/(2r)`. The correct BIE substitution is `D_s+S(q_s/Phi)=0`; production code uses the opposite S sign and omits `1/(2r)`.
8. For `R(z)=R_b+R'(z+h)`, inward `n_s=(-e_r+R'e_z)/sqrt(1+R'^2)`, so
   `q_s/Phi=[i kappa+1/(2r)+R' kappa tanh(kappa(z+h))]/sqrt(1+R'^2)`.
   This formula is correct only when `kappa` is the horizontal (generally complex) root consistent with the local FS Robin coefficient and bottom condition. In a sponge plateau, substituting the undamped real physical `k` is inconsistent.
9. Linear Bernoulli pressure is `p=-iρ ω Phi`. Force on the body is `Q=-integral(p n_B dS)=+iρ ω integral(Phi n_B dS)`. Production radiation and excitation integrations use the opposite body-force sign. The radiation RHS and force errors can double-cancel in A/B, so preserving today's positive-looking coefficients is not evidence either is correct.
10. With physical `Q_R=(ω²A-iωB)xi`, extraction is `A=Re(Q_R)/(ω²xi)` and `B=-Im(Q_R)/(ω xi)`. The RAO matrix `-ω²(M+A)+iωB+C` is correct for `exp(+iωt)` **if** A, B and excitation are physical body-force quantities. Current excitation sign gives a 180-degree phase error even where amplitude is unchanged.
11. Direct flux must use a pre-sponge outward control normal. For diagonal/unit displacement,
    `P_out=1/2 Re integral(p q_out* dS)=1/2 ω²B`, hence
    `B=(ρ/ω) Im integral(Phi q_out* dS)`.
    For cross modes, build the correctly indexed Hermitian pressure-work form and retain its raw non-Hermitian/non-real residual before reporting a real symmetric B. Do not PSD-clip it.
12. Existing “energy” is Haskind-derived, not flux. The current Haskind plus-sign integrand is not the antisymmetric Green concomitant under one normal convention and remains `NOT VALIDATED` pending a manufactured identity.

Mandatory minimal-test/fix order:

1. **T1 time/amplitude** — incident phase propagation, dynamic/kinematic FS identities, peak amplitude and pressure. This must pass without solver fitting.
2. **T2 Green/source orientation** — constant harmonic field on a simple watertight closed boundary using inward source normals; verify `D_s 1=0`, normal reversal, winding reversal and analytic off-surface derivative. This gates all matrix changes.
3. **T3/T4 physical FS and bottom** — finite-depth plane mode must satisfy `partial_z Phi=(ω²/g)Phi` at z=0 and zero bottom derivative; separately test the matrix substitution signs.
4. **T5 cylinder first, then frustum** — manufactured Hankel/finite-depth mode tests `ikappa`, `1/(2r)`, BIE S sign and the vertical projection. Include a complex-kappa sponge case derived from its Robin coefficient. This gates factors 2/4 and all outer studies.
5. **T6 unbounded translating sphere radiation** — exact dipole potential checks body `q_s=iωn_j` and the negative radiation RHS without FS/outer complications.
6. **T10 pressure traction, then T7 coefficient extractor** — analytic panel traction establishes `Q=-integral p n_B`; exact translating-sphere added mass verifies the coordinated RHS/force/A/B convention. Check before/after signs, not magnitudes only.
7. **T8 sphere diffraction** — manufactured/analytic normal cancellation verifies `q_D=-q_I` and excitation body-force sign.
8. **T9 direct flux** — a finite-depth progressive plane wave must reproduce `E c_g` with peak-amplitude factor 1/2; then test radius/quadrature convergence and cross-mode raw/Hermitian bookkeeping.
9. **T11 Haskind last** — only after pressure and flux pass; use a manufactured Green identity to determine normal, sign, conjugation, incident direction and prefactor. Until then report `B_haskind: NOT VALIDATED`.

Production experiments and the frustum factors remain blocked until tests 1–6 pass and the corresponding convention defects are fixed one coherent convention issue at a time. Flux/Haskind conclusions additionally require tests 8/9 in the list above.

### Checkpoint 4 — Production-code modification

Review the focused diff plus pre-change baseline, theory, minimal reproduction/test, clean representative reruns, convergence impact, and absence of forbidden shortcuts. Insufficient evidence is a `BLOCKER`.

**Current verdict:** `PENDING — first atomic fix authorized; no post-change diff/evidence reviewed`.

#### CP4 pre-change authorization 1 — double-layer source orientation

**Decision:** `PASS — authorization limited to this one issue`.

Accepted fail-before evidence:

- All 19 rows match their declared before-state, while physical `currentPass` remains visibly 10/19; expected failures are not counted as validation.
- T2 uses the actual `rankine_panel_integrals` kernel. The all-inward closed cube satisfies `D_s 1=0` to `5.55e-17`; changing only the top/FS to current outward orientation produces residual `0.21795`. Both winding-alignment controls pass, isolating a mixed domain/source orientation defect.
- T3/T5/T6/T10/T7 independently remain red with the expected signs/magnitudes. The complex-kappa row is honestly `NOT_IMPLEMENTED`.

Authorized implementation contract:

1. Canonical double-layer/source orientation is `n_D=-n_Omega` for every component. Mapping is body `+1`, FS `-1`, bottom `+1`, outer `+1` relative to the stored/winding derivative.
2. Apply the component sign to each source-column `dGdn` contribution (after image summation is equivalent) before forming `Dij=-dGdn/(4pi)`. It must affect **all collocation rows** for that source column.
3. Do not multiply `Sij`, prescribed Neumann RHS, stored physical normals, pressure normals, or the `+0.5` boundary self coefficient. The self coefficient is added once and is orientation-independent after canonicalization.
4. Assert stored normal–winding alignment per component before using the fixed component map. A silent inconsistent BMF must fail rather than receive the wrong canonical sign.
5. Store the named convention and component signs in `assembly_info` and the effective/cache manifest. The `assemble_rankine_matrix.m` file hash already invalidates potential cache through `get_rankine_code_version`; any new shared production helper must also be added to that version fingerprint.
6. Preserve fail-before CSV/MAT/MD unchanged and create separate after evidence. T2 mixed case must turn green near machine precision; T3/T5/T6/T10/T7 must remain in their prior physical fail state, showing this patch did not smuggle other convention fixes.
7. The after-test must touch the production canonicalization path. A hard-coded second copy of the orientation multiplication is not sufficient by itself. Diff inspection plus production-reported signs/hash are mandatory.
8. No WAMIT comparison or outer-domain solve is evidence for this patch. No reciprocity projection, PSD clipping, scaling, smoothing or unrelated modification is authorized.

### Checkpoint 5 — Final acceptance

Reconcile every checklist item against actual artifacts/data, answer Q1–Q6, audit conclusions against raw data, enumerate unresolved blockers, and emit final Auditor state variables. Only `FINAL_AUDIT = PASS` permits Phase-complete language.

**Current verdict:** `PENDING`.

## 11. Auditor risk watchlist

1. Stale result/potential cache or a cache key that does not include the varied factor.
2. Config changed but generated active panels/BC coefficients did not; mesh hashes or matrix rows stay unexpectedly identical.
3. A “resolution” run also moves radii, sponge start, zone width, or other component spacings.
4. Nominal spacing reported instead of actual mesh statistics and `h/lambda` at each frequency.
5. Small linear algebra residual mistaken for accurate boundary physics; collocation residual circularly restates the solved rows.
6. Pressure/flux/Haskind paths share arrays/prefactors or one is overwritten, so they are not independent.
7. Reciprocity projection, PSD clipping, smoothing, or selected plots hide raw signs/eigenvalues/outliers.
8. Closeness to WAMIT is labelled convergence or used to choose numerical parameters.
9. Condition estimates are unreliable/unlabelled, or worsening refinement is rationalized without local residual/potential evidence.
10. Final causal language exceeds experiment identifiability; an unresolved mixed cause is presented as a unique root cause.

## 12. Decision log

| Time (Asia/Shanghai) | Checkpoint | Decision | Summary |
|---|---|---|---|
| 2026-08-30 | Checklist initialization | PASS | Requirements captured; latest numerical audit read; production tree intentionally untouched by Auditor. Awaiting Executor plan for Checkpoint 1. |
| 2026-08-30 | Checkpoint 1 — Executor plan | PASS | Plan covers fixed Fine-body clean baseline at omega {0.5, 1.1, 1.5, 1.7, 2.0}, all nine OFAT factors, effective-geometry/cache provenance, wavelength normalization, three radius×resolution groups, `exp(+iωt)` BC audit, three independent damping paths, raw/local diagnostics, conditioning, WAMIT-only validation, single-change evidence gates, and required deliverables. Six carry-forward warnings are explicit under Checkpoint 1; none currently blocks design work. |
| 2026-08-30 | Checkpoint 2 — Outer-domain design | BLOCKER | Clean/cache authenticity, directional `h/lambda`, sponge start/width independence, and radius×resolution groups pass design review. Replacing physical FS/bottom radial extents with fine-resolution-zone endpoints does not meet two explicit OFAT requirements. A valid independently controlled closed-domain parameterization plus derived connector BC, or an explicit unresolved requirement/scope decision, is required. |
| 2026-08-30 | Checkpoint 2 resubmission — Outer-domain design | PASS | The physical extent gap is resolved by independent top/bottom frustum radii with deterministic closure effects. Derived meridional/tangential counts preserve actual wall spacing and prevent a radius-resolution confound. Frustum topology/metric assertions are explicit. Its projected finite-depth outgoing BC remains a mandatory Checkpoint 3 gate before factors 2/4 may run. |
| 2026-08-30 | Checkpoint 3 — BC mathematics | BLOCKER | Independent Green derivation confirms `D_s=1/2-K_s` requires inward source normals and RHS `-S q_s`. Current FS orientation, radiation RHS, finite-depth FS coefficient, outer S sign/decay/sponge consistency, and body-force signs are inconsistent; radiation RHS and force errors can double-cancel. Diffraction body RHS, bottom zero flux, incident time convention and RAO matrix form are analytically consistent. Prioritized T1–T11 gates are recorded; no production change or outer experiment is yet accepted. |
| 2026-08-30 | CP4 pre-change authorization 1 | PASS | Fail-before suite truthfully reports 10/19 physical passes and isolates mixed FS double-layer orientation with actual Rankine kernel/closed-cube controls. Authorized only per-source-column FS `dGdn` canonicalization plus alignment assertions and cache/assembly metadata. Separate production-path after evidence must turn only T2 green; full CP4 remains pending. |
