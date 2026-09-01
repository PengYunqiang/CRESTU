function [force, pressure] = compute_wave_excitation(phi_incident, phi_diffraction, nj, areas, omega, rho)
% COMPUTE_WAVE_EXCITATION Integrate total first-order pressure into generalized wave-excitation loads.
%
% Syntax:
%   [force, pressure] = compute_wave_excitation(phi_incident, phi_diffraction, nj, areas, omega, rho)
%
% Description:
%   Computes hydrodynamic coefficients, loads, or rigid-body response.
%   Results retain the CRESTU global 6-DOF order.
%
% Inputs:
%   phi_incident       - [N x Nh] Complex incident-wave potentials, [m^2/s].
%   phi_diffraction    - [N x Nh] Complex diffraction potentials, [m^2/s].
%   nj                 - [N x 6Nb] Generalized panel normals, with rotational columns in [m].
%   areas              - [N x 1] Panel areas, [m^2].
%   omega              - [scalar] Angular frequency, [rad/s].
%   rho                - [scalar] Fluid density, [kg/m^3].
%
% Outputs:
%   force              - [Ndof x Nh] Complex generalized hydrodynamic loads, [N] and [N m].
%   pressure           - [N x Nh] Complex first-order dynamic pressure, [Pa].
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
    if isempty(phi_incident)
        phi_incident = zeros(size(phi_diffraction),'like', phi_diffraction);
    end
    if isempty(phi_diffraction)
        phi_diffraction = zeros(size(phi_incident),'like', phi_incident);
    end
    if ~isequal(size(phi_incident), size(phi_diffraction)) || ...
            size(phi_incident, 1) ~= size(nj, 1) || numel(areas) ~= size(nj, 1)
        error('CRESTU:ExcitationShape','Incident/diffraction potentials, normals, and areas are inconsistent.');
    end
    pressure = -1i * omega * rho * (phi_incident + phi_diffraction);
% Body generalized load is the negative pressure traction integral.
    force = -nj.' * (pressure .* areas);
end
