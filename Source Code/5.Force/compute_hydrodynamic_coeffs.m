function [added_mass,damping,diagnostics] = compute_hydrodynamic_coeffs(phi_radiation,nj,areas,omega,rho,varargin)
% COMPUTE_HYDRODYNAMIC_COEFFS Execute the documented compute_hydrodynamic_coeffs operation.
%
% Syntax:
%   [added_mass,damping,diagnostics] = compute_hydrodynamic_coeffs(phi_radiation,nj,areas,omega,rho,varargin)
%
% Inputs:
%   phi_radiation   : [N x Ndof] Complex radiation potential, in m^2/s.
%   nj              : [N x 6Nb] Generalized panel normals for all body modes.
%   areas           : [N x 1] Panel areas, in m^2.
%   omega           : [scalar] Angular frequency, in rad/s.
%   rho             : [scalar] Water density, in kg/m^3.
%   varargin        : [documented value] Input required by the implemented function contract.
%
% Outputs:
%   added_mass      : [Ndof x Ndof x Nf] Added-mass matrices in SI translational/rotational units.
%   damping         : [Ndof x Ndof x Nf] Radiation-damping matrices in SI translational/rotational units.
%   diagnostics     : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
%COMPUTE_HYDRODYNAMIC_COEFFS Added mass and damping for exp(i*omega*t).
% Radiation BC: d(phi_j)/dn=i*omega*n_j. Fluid force is
% H_ij=int(-i*omega*rho*phi_j)n_i dS = omega^2*A_ij-i*omega*B_ij.
    validateattributes(omega,{'numeric'},{'scalar','real','positive','finite'});
    validateattributes(rho,{'numeric'},{'scalar','real','positive','finite'});
    areas=reshape(areas,[],1);
    if size(phi_radiation,1)~=size(nj,1)||numel(areas)~=size(nj,1)
        error('CRESTU:ForceShape','Potential, generalized normals, and areas have inconsistent rows.');
    end
    if size(phi_radiation,2)~=size(nj,2)
        error('CRESTU:RadiationModes','Radiation potential must have one column per generalized mode.');
    end
    potential_integrals=nj.'*(phi_radiation.*areas);
    symmetry_weights=ones(size(potential_integrals));
    if ~isempty(varargin)
        symmetry=varargin{1};
        if isstruct(symmetry)&&isfield(symmetry,'mode_parity')
            symmetry_weights=symmetry_force_weights(symmetry.mode_parity, ...
                symmetry.mode_parity,symmetry.isx,symmetry.isy);
            potential_integrals=potential_integrals.*symmetry_weights;
        end
    end
    radiation_force=-1i*omega*rho*potential_integrals;
    added_mass_raw=real(radiation_force)/(omega^2);
    damping_raw=-imag(radiation_force)/omega;
    scaleA=max(norm(added_mass_raw,'fro'),eps); scaleB=max(norm(damping_raw,'fro'),eps);
    raw_added_mass_symmetry_error=norm(added_mass_raw-added_mass_raw.','fro')/scaleA;
    raw_damping_symmetry_error=norm(damping_raw-damping_raw.','fro')/scaleB;

    % Enforce the reciprocal projection after retaining the raw residual as
    % a mesh/solver diagnostic. This is the nearest Frobenius-norm matrix
    % satisfying the potential-flow reciprocity identities A=A' and B=B'.
    added_mass=0.5*(added_mass_raw+added_mass_raw.');
    damping=0.5*(damping_raw+damping_raw.');
    diagnostics=struct('radiation_force',radiation_force, ...
        'raw_added_mass_symmetry_error',raw_added_mass_symmetry_error, ...
        'raw_damping_symmetry_error',raw_damping_symmetry_error, ...
        'added_mass_symmetry_error',norm(added_mass-added_mass.','fro')/scaleA, ...
        'damping_symmetry_error',norm(damping-damping.','fro')/scaleB, ...
        'min_added_mass_diagonal',min(real(diag(added_mass))), ...
        'min_damping_diagonal',min(real(diag(damping))), ...
        'symmetry_weights',symmetry_weights);
end
