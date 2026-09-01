function [added_mass, damping, diagnostics] = compute_hydrodynamic_coeffs(phi_radiation, nj, areas, omega, rho, varargin)
% COMPUTE_HYDRODYNAMIC_COEFFS Extract added mass and radiation damping from radiation potentials.
%
% Syntax:
%   [added_mass, damping, diagnostics] = compute_hydrodynamic_coeffs(phi_radiation, nj, areas, omega, rho, varargin)
%
% Description:
%   Computes hydrodynamic coefficients, loads, or rigid-body response.
%   Results retain the CRESTU global 6-DOF order.
%
% Inputs:
%   phi_radiation      - [N x Ndof] Complex radiation potentials, [m^2/s].
%   nj                 - [N x 6Nb] Generalized panel normals, with rotational columns in [m].
%   areas              - [N x 1] Panel areas, [m^2].
%   omega              - [scalar] Angular frequency, [rad/s].
%   rho                - [scalar] Fluid density, [kg/m^3].
%   varargin           - [cell array] Optional algorithm controls documented by the function implementation.
%
% Outputs:
%   added_mass         - [Ndof x Ndof] Added-mass matrix in consistent translational and rotational SI units.
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

%% Stage 1: Validate Inputs and Initialize the Algorithm

    validateattributes(omega, {'numeric'}, {'scalar','real','positive','finite'});
    validateattributes(rho, {'numeric'}, {'scalar','real','positive','finite'});
    areas = reshape(areas, [], 1);
    if size(phi_radiation, 1) ~= size(nj, 1) || numel(areas) ~= size(nj, 1)
        error('CRESTU:ForceShape','Potential, generalized normals, and areas have inconsistent rows.');
    end
    if size(phi_radiation, 2) ~= size(nj, 2)
        error('CRESTU:RadiationModes','Radiation potential must have one column per generalized mode.');
    end
    % <<<CORE>>> hydrodynamic_coefficient_projection, paper_eq=none, benchmark=single_sphere_reciprocity
    potential_integrals = nj.' * (phi_radiation .* areas);
    symmetry_weights = ones(size(potential_integrals));
    if ~isempty(varargin)
        symmetry = varargin{1};
        if isstruct(symmetry) && isfield(symmetry,'mode_parity')
            symmetry_weights = symmetry_force_weights(symmetry.mode_parity, ...
                symmetry.mode_parity, symmetry.isx, symmetry.isy);
            potential_integrals = potential_integrals .* symmetry_weights;
        end
    end
% Body generalized load uses traction -p*n_B, with p=-i*rho*omega*phi.
    radiation_force = 1i * omega * rho * potential_integrals;
    added_mass_raw = real(radiation_force) / (omega^2);
    damping_raw = -imag(radiation_force) / omega;
    scaleA = max(norm(added_mass_raw,'fro'), eps);
    scaleB = max(norm(damping_raw,'fro'), eps);
    raw_added_mass_symmetry_error = norm(added_mass_raw - added_mass_raw.','fro') / scaleA;
    raw_damping_symmetry_error = norm(damping_raw - damping_raw.','fro') / scaleB;

% Enforce the reciprocal projection after retaining the raw residual as
% a mesh/solver diagnostic. This is the nearest Frobenius-norm matrix
% satisfying the potential-flow reciprocity identities A=A' and B=B'.
    added_mass = 0.5 * (added_mass_raw + added_mass_raw.');
    damping = 0.5 * (damping_raw + damping_raw.');
    symmetric_damping_raw = 0.5 * (damping_raw + damping_raw.');
    minimum_eigenvalue_symmetric_damping_raw = ...
        min(real(eig(symmetric_damping_raw)));
    diagnostics = struct('radiation_force', radiation_force, ...
'added_mass_raw', added_mass_raw, 'damping_raw', damping_raw, ...
'added_mass_reciprocal', added_mass, 'damping_reciprocal', damping, ...
'raw_added_mass_symmetry_error', raw_added_mass_symmetry_error, ...
'raw_damping_symmetry_error', raw_damping_symmetry_error, ...
'added_mass_symmetry_error', norm(added_mass - added_mass.','fro') / scaleA, ...
'damping_symmetry_error', norm(damping - damping.','fro') / scaleB, ...
'min_added_mass_diagonal', min(real(diag(added_mass))), ...
'min_damping_diagonal', min(real(diag(damping))), ...
'added_mass_raw_diagonal', diag(added_mass_raw), ...
'damping_raw_diagonal', diag(damping_raw), ...
'symmetric_damping_raw', symmetric_damping_raw, ...
'minimum_eigenvalue_symmetric_damping_raw', ...
minimum_eigenvalue_symmetric_damping_raw, ...
'symmetry_weights', symmetry_weights);
    % <<</CORE>>>
end
