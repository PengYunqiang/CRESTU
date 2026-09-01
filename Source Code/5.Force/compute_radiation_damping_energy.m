function [damping, diagnostics] = compute_radiation_damping_energy( ...
        phi_radiation, nj, centers, normals, areas, omega, cfg, mode_parity, n_theta)
% COMPUTE_RADIATION_DAMPING_ENERGY Deprecated compatibility wrapper only.
% Active scientific analysis must call compute_haskind_damping_candidate.
%
% Syntax:
%   [damping, diagnostics] = compute_radiation_damping_energy(phi_radiation, nj, centers, normals, areas, omega, cfg, mode_parity, n_theta)
%
% Description:
%   Computes hydrodynamic coefficients, loads, or rigid-body response.
%   Results retain the CRESTU global 6-DOF order.
%
% Inputs:
%   phi_radiation      - [N x Ndof] Complex radiation potentials, [m^2/s].
%   nj                 - [N x 6Nb] Generalized panel normals, with rotational columns in [m].
%   centers            - [N x 3] Panel collocation points in global coordinates, [m].
%   normals            - [N x 3] Unit panel normals, dimensionless.
%   areas              - [N x 1] Panel areas, [m^2].
%   omega              - [scalar] Angular frequency, [rad/s].
%   cfg                - [struct] Validated CRESTU configuration containing SI-valued physical and numerical parameters.
%   mode_parity        - [Ndof x 2] Reflection parity signs for the rigid-body modes, dimensionless.
%   n_theta            - [scalar] Number of azimuthal integration samples, dimensionless.
%
% Outputs:
%   damping            - [Ndof x Ndof] Radiation-damping matrix in consistent translational and rotational SI units.
%   diagnostics        - [struct] Numerical conditioning, symmetry, residual, or reconstruction diagnostics.
%
% Governing Equations / Theory:
%   Linear unsteady Bernoulli pressure, generalized surface integration, hydrostatics, radiation energy, or the frequency-domain rigid-body equation as applicable.
%
% References:
%   - Newman, J. N. (1977), Marine Hydrodynamics; Faltinsen, O. M. (1990), Sea Loads on Ships and Offshore Structures.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

    if nargin < 9
        n_theta = [];
    end
    [damping, candidateDiagnostics] = ...
        compute_haskind_damping_candidate(phi_radiation, nj, centers, ...
        normals, areas, omega, cfg, mode_parity, n_theta);
    diagnostics = struct('deprecatedCompatibility', struct( ...
        'legacyFunctionName', 'compute_radiation_damping_energy', ...
        'replacementFunctionName', ...
        'compute_haskind_damping_candidate', ...
        'semanticWarning', ['The legacy name is not an energy-flux path; ', ...
        'consume candidateDiagnostics only as NOT_VALIDATED.'], ...
        'candidateDiagnostics', candidateDiagnostics));
end
