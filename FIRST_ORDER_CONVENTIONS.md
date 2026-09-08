# First-order conventions and corrected numerical path

This document is normative for the supplied active source. Older `TheoryManual.md` prose and historical reports contain incompatible force/Haskind signs. They are evidence of past versions, not authority over the following derivation. This revision has Python auxiliary tests; it is **not yet MATLAB runtime validation**.

## Harmonic, normal and potential definitions

All dimensional fields are real parts of complex amplitudes times `exp(+i*omega*t)`. Velocity is grad(Phi), dynamic pressure is `p = -i*rho*omega*Phi`. Body normals `n_B` point OUT of the body, INTO the fluid. The generalized normals are `N=[n_B, (x-CG) cross n_B]`, in surge, sway, heave, roll, pitch, yaw order per body. A positive generalized force is `F=-integral(p*N)dS`.

`Phi_R,j` is the potential per unit generalized displacement (metre for translation, radian for rotation); hence its body derivative is `q_B,j=i*omega*N_j`. The unit-velocity potential is `psi_j=Phi_R,j/(i*omega)`. These two definitions must not be interchanged without the frequency factor.

WAMIT's theoretical normal points OUT of the fluid and is therefore `n_W=-n_B` on the body. Its unit-velocity potential solves the same physical problem when both derivative and generalized normal are transformed. An isolated minus sign from the WAMIT force formula cannot be copied into a body-outward implementation.

## Rankine direct boundary integral equation

Use `G=1/|x-y|` and divide its panel integrals by `4*pi`. Convert the stored normals to one canonical inward-fluid normal `n_s`: body +1, free surface -1, bottom +1, outer wall +1 for this mesh convention. Green's representation is

`Phi(x) = integral[(dG/dn_s)*Phi - G*q_s]dS/(4*pi)`.

At a smooth planar collocation panel,

`(0.5*I - D_s)*Phi = -S*q_s`.

Thus radiation RHS is `-S_body*(i*omega*N)`; diffraction/scattering RHS is `+S_body*q_I`, because `q_scattered=-q_I`. No extra 2, pi, rho or omega belongs in either RHS. The self double-layer principal value is zero; the solid-angle jump is the explicit diagonal 0.5.

Warped quads are approximated by the plane through their centroid with diagonal-cross normal (the mid-edge plane). BOTH single- and double-layer formulas must use the vertices projected to that plane. The original source projected only the edge part but used nonplanar vertices for the solid angle. The patch restores `dG/dn_source = -dG/dn_field`. Small inter-panel gaps are an acknowledged low-order planar approximation, to be reduced by body refinement.

## Surface, bottom and outer boundaries

Without absorption, `dPhi/dz=(omega^2/g)*Phi` at z=0, regardless of finite depth. The real propagating wavenumber solves `k*tanh(k*h)=omega^2/g`; replacing the FS coefficient by k in finite depth is wrong. In the current quadratic sponge, `nu=(omega^2/g)*(1-i*mu(r))`, so canonical `q_s=-nu*Phi` and the matrix column is `D_solved - nu*S`.

For `exp(+i*omega*t)`, an outgoing cylindrical wave is proportional to `H_m^(2)(k*r)`, asymptotically `exp(-i*k*r)/sqrt(r)`. Thus an inward radial normal gives `q_s≈(i*k+1/(2*r))*Phi`. The current complex k branch has positive real and nonpositive imaginary parts and satisfies the damped finite-depth dispersion relation. The finite wall and sponge remain approximations requiring domain convergence; neither is a justification to change production load signs.

## Incident, diffraction, loads and response

For unit surface elevation at the origin and heading beta,

`Phi_I = i*g/omega * cosh(k*(z+h))/cosh(k*h) * exp(-i*k*(x*cos(beta)+y*sin(beta)))`.

`Phi_total = Phi_I + Phi_scattered`. Total excitation is

`F = +i*rho*omega*integral(N*Phi_total)dS`.

Define `I_ij=integral(N_i*Phi_R,j)dS`. Then

`F_R,ij=+i*rho*omega*I_ij = omega^2*A_ij - i*omega*B_ij`;

`A=real(F_R)/omega^2`, `B=-imag(F_R)/omega`.

The displacement equation is `[-omega^2*(M+A)+i*omega*B+C]*xi=F`. Mass, hydrostatics and pressure damping are not adjusted to fit RAO.

## Haskind and energy are diagnostics, not replacement damping

Reciprocity between unit-velocity radiation and scattering gives

`integral(N_j*Phi_scattered)dS = -integral(psi_j*q_I)dS`.

Therefore, in this normal convention,

`F_H,j = +i*rho*omega*integral(N_j*Phi_I - psi_j*q_I)dS`.

The previous Haskind helper had the opposite exterior sign and a plus instead of the internal minus. This is a genuine bug in an auxiliary path. Correcting it does not authorize substituting its damping for production B.

For outgoing radiation per unit displacement on a control cylinder in the undamped region,

`B_flux,jj = rho/omega * integral(Im(Phi_R,j*conj(dPhi_R,j/dr)))dS`.

The Haskind angular-energy matrix has prefactor `k/(8*pi*rho*g*c_g)` multiplying `integral(F_H(beta)*F_H(beta)^H)d beta`, with `c_g=domega/dk`. Check its real symmetric part against pressure damping and control-cylinder flux. A finite discretized absorbing-domain solve need not satisfy these identities exactly; examine convergence rather than swapping formulas.

## Reference normalization and provenance

The supplied FullSphereIRR3 GDF has ULEN=1 m, diameter 10 m, water depth 50 m, rho=1025 kg/m3 and g=9.80665 m/s2. Its 300 input panels are coordinate-identical to the supplied Medium body, not the 588-panel Fine body. The supplied deck requests angular frequencies, but `.1/.2/.4` column 1 is OUTPUT PERIOD IN SECONDS. Match `2*pi/period` independently for each extension.

Let r_i=0 for each body's local translation DOF and 1 for each local rotation DOF. Then `A=rho*ULEN^(3+r_i+r_j)*Ahat`, `B=rho*omega*ULEN^(3+r_i+r_j)*Bhat`, `F=rho*g*a*ULEN^(2+r_i)*Fhat`; per-unit-wave-amplitude RAO is `xihat/ULEN^r_i`. Option `.2` is Haskind TOTAL excitation; option `.3` is diffraction-pressure TOTAL excitation, not the isolated scattered component.

## Discretization and reproducibility

The original BodyMesh Fine had only 3 near-field and 2 sponge radial FS layers, and only the corresponding coarse bottom rings, despite automatic radius expansion to at least 1.5 wavelengths of sponge. A body refinement is not a free-surface or bottom refinement. The corrected supplied C/M/F study uses each body's own compatible waterline (24/40/56 sectors), FS 24+64 layers, bottom 6+12 layers and 2 vertical wall layers. It is therefore a coupled body/azimuthal refinement, not a body-only OFAT. The old fixed Fine waterline was incompatible with the coarser body polygons and left a seam gap; the new compatibility check rejects that input. These are explicit tested benchmark discretizations, not universal guaranteed settings. `audit_rankine_wave_resolution` reports actual sizes and warns without stopping the investigation. Refine again for other geometry/frequency ranges.

All source dependencies are now fingerprinted. Never validate current code from a historical MAT/CSV or manually copied old cache. The reproduction wrapper's reload mode reloads only a current clean-run cache. See `OUTPUT_FROM_GPT6PRO` for evidence and pending MATLAB tests.

## Primary sources

WAMIT official User Manual, Chapter 3 (definitions and normalization), Chapter 4 (input/output options), Chapter 15.1–15.3 (BVP, normals, Green identity and planar projection):

https://www.wamit.com/manual7.x/html/wamit_v75manualch3.html

https://www.wamit.com/manual7.x/html/wamit_v75manualch4.html

https://www.wamit.com/manualv7.4/wamit_v74manualch15.html
