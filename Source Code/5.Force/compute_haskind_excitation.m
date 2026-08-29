function force = compute_haskind_excitation(phi_radiation, nj, centers, normals, areas, omega, cfg, headings, mode_parity)
% COMPUTE_HASKIND_EXCITATION Evaluate wave-excitation loads from radiation solutions through Haskind reciprocity.
%
% Syntax:
%   force = compute_haskind_excitation(phi_radiation, nj, centers, normals, areas, omega, cfg, headings, mode_parity)
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
%   headings           - [1 x Nh] Wave headings, [deg].
%   mode_parity        - [Ndof x 2] Reflection parity signs for the rigid-body modes, dimensionless.
%
% Outputs:
%   force              - [Ndof x Nh] Complex generalized hydrodynamic loads, [N] and [N m].
%
% Governing Equations / Theory:
%   Linear unsteady Bernoulli pressure, generalized surface integration, hydrostatics, radiation energy, or the frequency-domain rigid-body equation as applicable.
%
% References:
%   - Newman, J. N. (1977), Marine Hydrodynamics; Faltinsen, O. M. (1990), Sea Loads on Ships and Offshore Structures.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    areas = reshape(areas, [], 1);
    headings = reshape(headings, 1, []);
    ndof = size(nj, 2);
    nh = numel(headings);
    force = complex(zeros(ndof, nh));
    if nargin < 9 || isempty(mode_parity)
        mode_parity = get_mode_parities(cfg.n_bodies, cfg.isx, cfg.isy);
    end
    psi = phi_radiation / (1i * omega);
    for j = 1:ndof
        parity = mode_parity(j, :);
        for headingIndex = 1:nh
            [phi_I, dphi_I] = decompose_incident_wave_symmetry(centers, normals, omega, cfg.grav, ...
                cfg.water_depth, headings(headingIndex), 1, cfg.isx, cfg.isy, parity);
            integral = sum((nj(:, j) .* phi_I + psi(:, j) .* dphi_I) .* areas);
            force(j, headingIndex) = -1i * omega * cfg.rho * cfg.symmetry.multiplicity * integral;
        end
    end
end
