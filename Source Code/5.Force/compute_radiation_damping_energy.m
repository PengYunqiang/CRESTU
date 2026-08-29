function [damping, diagnostics] = compute_radiation_damping_energy( ...
        phi_radiation, nj, centers, normals, areas, omega, cfg, mode_parity, n_theta)
% COMPUTE_RADIATION_DAMPING_ENERGY Evaluate radiation damping from the far-boundary energy flux.
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

%% Stage 1: Validate Inputs and Initialize the Algorithm

    if nargin < 9 || isempty(n_theta)
        n_theta = 72;
    end
    % <<<CORE>>> radiation_damping_energy_check, paper_eq=haskind_energy_identity, benchmark=single_sphere_pressure_energy
    headings = (0:n_theta - 1) * (360 / n_theta);
    dtheta = 2 * pi / n_theta;
    force = compute_haskind_excitation(phi_radiation, nj, centers, normals, areas, ...
        omega, cfg, headings, mode_parity);
    [wave_number, ~] = solve_wave_dispersion(omega, cfg.grav, cfg.water_depth);
    if cfg.water_depth > 0
        kh = wave_number * cfg.water_depth;
        depth_factor = tanh(kh) + kh / (cosh(kh)^2);
    else
        depth_factor = 1;
    end
    prefactor = omega * wave_number / (4 * pi * cfg.rho * cfg.grav^2 * depth_factor);
    damping = real(prefactor * (force * force') * dtheta);
    damping = 0.5 * (damping + damping.');
    eigenvalues = eig(damping);
    tolerance = 1e-10 * max(1, max(abs(eigenvalues)));
    if min(eigenvalues) < -tolerance
        warning('CRESTU:EnergyDampingPSD','Energy damping has minimum eigenvalue %g.', min(eigenvalues));
    end
    diagnostics = struct('method','Haskind energy flux','headings', headings, ...
'haskind_excitation', force,'depth_factor', depth_factor,'wavenumber', wave_number, ...
'min_eigenvalue', min(eigenvalues));
    % <<</CORE>>>
end
