function [phi_I, dphi_I_dn, k] = compute_incident_wave(centers, normals, omega, g, depth, beta_deg, amplitude)
% COMPUTE_INCIDENT_WAVE Evaluate the finite- or infinite-depth incident-wave potential and normal velocity.
%
% Syntax:
%   [phi_I, dphi_I_dn, k] = compute_incident_wave(centers, normals, omega, g, depth, beta_deg, amplitude)
%
% Description:
%   The routine implements a component of the linear Rankine boundary-element formulation for incompressible, irrotational gravity-wave flow. Geometry, reflection parity, free-surface impedance, and complex phase follow the project convention exp(i*omega*t).
%
% Inputs:
%   centers            - [N x 3] Panel collocation points in global coordinates, [m].
%   normals            - [N x 3] Unit panel normals, dimensionless.
%   omega              - [scalar] Angular frequency, [rad/s].
%   g                  - [scalar] Gravitational acceleration, [m/s^2].
%   depth              - [scalar] Water depth; a nonpositive value denotes infinite depth, [m].
%   beta_deg           - [scalar or vector] Wave-propagation heading, [deg].
%   amplitude          - [scalar] Incident-wave amplitude, [m].
%
% Outputs:
%   phi_I              - [N x Nh] Complex incident-wave potential, [m^2/s].
%   dphi_I_dn          - [N x Nh] Complex incident-wave normal velocity, [m/s].
%   k                  - [scalar] Gravity-wave wavenumber, [1/m].
%
% Governing Equations / Theory:
%   Green third identity, the Rankine kernel 1/r, linearized free-surface theory, and reflection symmetry as applicable.
%
% References:
%   - Newman, J. N. (1977), Marine Hydrodynamics; Hess and Smith (1964); project boundary-condition specification.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% --- 1. Validate Inputs and Initialize the Algorithm ---

    if nargin < 7 || isempty(amplitude), amplitude = 1; end
    centers = reshape(centers, [], 3); normals = reshape(normals, [], 3);
    if size(centers, 1) ~= size(normals, 1), error('CRESTU:IncidentShape', 'Centers/normals row mismatch.'); end
    validateattributes(omega, {'numeric'}, {'scalar', 'real', 'positive', 'finite'});
    k = solve_dispersion(omega, g, depth); beta = beta_deg * pi / 180;
    phase = exp(-1i * k * (centers(:, 1) * cos(beta) + centers(:, 2) * sin(beta)));
    coefficient = 1i * g * amplitude / omega;
    if depth > 0
        vertical = cosh(k * (centers(:, 3) + depth)) / cosh(k * depth);
        vertical_dz = k * sinh(k * (centers(:, 3) + depth)) / cosh(k * depth);
    else
        vertical = exp(k * centers(:, 3)); vertical_dz = k * vertical;
    end
    phi_I = coefficient * vertical .* phase;
    u = (-1i * k * cos(beta)) * phi_I; v = (-1i * k * sin(beta)) * phi_I;
    w = coefficient * vertical_dz .* phase;
    dphi_I_dn = u .* normals(:, 1) + v .* normals(:, 2) + w .* normals(:, 3);
end

function k = solve_dispersion(omega, g, depth)
% SOLVE_DISPERSION Solve the local linear gravity-wave dispersion relation.
%
% Syntax:
%   k = solve_dispersion(omega, g, depth)
%
% Description:
%   The routine implements a component of the linear Rankine boundary-element formulation for incompressible, irrotational gravity-wave flow. Geometry, reflection parity, free-surface impedance, and complex phase follow the project convention exp(i*omega*t).
%
% Inputs:
%   omega              - [scalar] Angular frequency, [rad/s].
%   g                  - [scalar] Gravitational acceleration, [m/s^2].
%   depth              - [scalar] Water depth; a nonpositive value denotes infinite depth, [m].
%
% Outputs:
%   k                  - [scalar] Gravity-wave wavenumber, [1/m].
%
% Governing Equations / Theory:
%   Green third identity, the Rankine kernel 1/r, linearized free-surface theory, and reflection symmetry as applicable.
%
% References:
%   - Newman, J. N. (1977), Marine Hydrodynamics; Hess and Smith (1964); project boundary-condition specification.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% --- 1. Validate Inputs and Initialize the Algorithm ---

    k = omega^2 / g;
    if depth <= 0, return; end
    k = max(k, omega / sqrt(g * depth));
    for iter = 1:50
        kh = k * depth; t = tanh(kh); f = g * k * t - omega^2;
        df = g * t + g * k * depth / (cosh(kh)^2); step = f / df; k_new = max(k - step, eps);
        if abs(k_new - k) <= 1e-12 * max(1, k), k = k_new; return; end
        k = k_new;
    end
    warning('CRESTU:DispersionConvergence', 'Dispersion iteration reached its limit at omega=%g.', omega);
end
