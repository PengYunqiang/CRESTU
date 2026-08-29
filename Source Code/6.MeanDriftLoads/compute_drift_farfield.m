function drift = compute_drift_farfield(mesh, state, cfg)
% COMPUTE_DRIFT_FARFIELD Evaluate Maruo-Newman mean drift loads from Kochin functions.
%
% Syntax:
%   drift = compute_drift_farfield(mesh, state, cfg)
%
% Description:
%   Computes or processes second-order mean wave-drift quantities.
%   Complex averages follow exp(i*omega*t) and the 6-DOF order.
%
% Inputs:
%   mesh               - [struct] Boundary-panel mesh with Cartesian geometry in SI units.
%   state              - [struct] Frequency-domain potential, motion, and load state in SI units.
%   cfg                - [struct] Validated CRESTU configuration containing SI-valued physical and numerical parameters.
%
% Outputs:
%   drift              - [struct] Mean second-order generalized loads and audited component terms, [N] and [N m].
%
% Governing Equations / Theory:
%   Pinkster near-field direct pressure integration, Maruo-Newman far-field momentum balance, or supporting surface reconstruction as applicable.
%
% References:
%   - Pinkster, J. A. (1980), Low Frequency Second Order Wave Exciting Forces on Floating Structures; Newman, J. N. (1974), second-order slowly varying forces.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    required = {'phi','dphi_dn','omega','headings'};
    for k = 1:numel(required)
        if ~isfield(state, required{k})
            error('CRESTU:DriftState','Far-field state lacks %s.', required{k});
        end
    end
    if cfg.isx || cfg.isy
        error('CRESTU:FarfieldSymmetryExpansion', ...
'Kochin drift currently requires full-domain body fields (ISX=ISY=0).');
    end
    np = mesh.n_panels;
    nh = size(state.phi, 2);
    areas = reshape(mesh.areas, np, 1);
    centers = reshape(mesh.centers, np, 3);
    normals = reshape(mesh.normals, np, 3);
    [k, ~] = solve_wave_dispersion(state.omega, cfg.grav, cfg.water_depth);
    nu = state.omega^2 / cfg.grav;
    if isfield(state,'n_theta')
        ntheta = state.n_theta;
    else
        ntheta = 360;
    end
    theta = linspace(0, 2 * pi, ntheta + 1).';
    theta(end) = [];
    dtheta = 2 * pi / ntheta;
    H = complex(zeros(ntheta, nh));
    for t = 1:ntheta
        ct = cos(theta(t));
        st = sin(theta(t));
        phase = exp(-1i * k * (centers(:, 1) * ct + centers(:, 2) * st));
        if cfg.water_depth > 0
            fz = cosh(k * (centers(:, 3) + cfg.water_depth)) / cosh(k * cfg.water_depth);
            dfz = k * sinh(k * (centers(:, 3) + cfg.water_depth)) / cosh(k * cfg.water_depth);
        else
            fz = exp(k * centers(:, 3));
            dfz = k * fz;
        end
        psi = fz .* phase;
        dpsi = psi .* (-1i * k * (normals(:, 1) * ct + normals(:, 2) * st)) + dfz .* phase .* normals(:, 3);
        for headingIndex = 1:nh
% CRESTU panel normals point from the body into the fluid. The
% exterior-fluid Green identity uses the opposite normal on the
% wetted-body boundary, hence the leading minus sign.
            H(t, headingIndex) = -sum((psi .* state.dphi_dn(:, headingIndex) - state.phi(:, headingIndex) .* dpsi) .* areas);
        end
    end
    Hp = (circshift(H, -1, 1) - circshift(H, 1, 1)) / (2 * dtheta);
    loads = zeros(6, nh);
    denominator = nu + cfg.water_depth * (k^2 - nu^2);
    if cfg.water_depth <= 0
        denominator = nu;
    end
    prefactor = cfg.rho * k^3 / (8 * pi * nu) * (nu / denominator);
% For a diagonal mean-drift pair, the two conjugate incident-scattered
% terms in OrcaWave equation (12) combine into twice the real part.
    cross_prefactor = cfg.rho * k * state.omega / (2 * nu);
    for headingIndex = 1:nh
        beta = state.headings(headingIndex) * pi / 180;
        loads(1, headingIndex) = real(prefactor * sum(cos(theta) .* H(:, headingIndex) .* conj(H(:, headingIndex))) * dtheta + ...
            cross_prefactor * cos(beta) * interp_periodic(theta, H(:, headingIndex), pi + beta));
        loads(2, headingIndex) = real(prefactor * sum(sin(theta) .* H(:, headingIndex) .* conj(H(:, headingIndex))) * dtheta + ...
            cross_prefactor * sin(beta) * interp_periodic(theta, H(:, headingIndex), pi + beta));
        yaw_integral = sum(Hp(:, headingIndex) .* conj(H(:, headingIndex))) * dtheta;
        yaw_cross = interp_periodic(theta, Hp(:, headingIndex), pi + beta);
        yaw_incident_pair = yaw_cross - conj(yaw_cross);
        loads(6, headingIndex) = real(-1i * cfg.rho * k^2 / (8 * pi * nu) * (nu / denominator) * yaw_integral + ...
            1i * cfg.rho * state.omega / (4 * nu) * yaw_incident_pair);
    end
    drift = struct('method','far_field_kochin','omega', state.omega,'headings', state.headings, ...
'theta', theta,'kochin', H,'kochin_derivative', Hp,'total', loads, ...
'limitations','Only global surge, sway and yaw; ensemble load for multi-body cases.');
end

function value = interp_periodic(theta, values, query)
% INTERP_PERIODIC Interpolate a sampled periodic angular function.
%
% Syntax:
%   value = interp_periodic(theta, values, query)
%
% Description:
%   Computes or processes second-order mean wave-drift quantities.
%   Complex averages follow exp(i*omega*t) and the 6-DOF order.
%
% Inputs:
%   theta              - [N x 1] Periodic angular coordinates, [rad].
%   values             - [numeric array] Samples to be transformed or interpolated; units are preserved.
%   query              - [scalar] Periodic interpolation coordinate, [rad].
%
% Outputs:
%   value              - [scalar] Interpolated or parsed value with units inherited from the input contract.
%
% Governing Equations / Theory:
%   Pinkster near-field direct pressure integration, Maruo-Newman far-field momentum balance, or supporting surface reconstruction as applicable.
%
% References:
%   - Pinkster, J. A. (1980), Low Frequency Second Order Wave Exciting Forces on Floating Structures; Newman, J. N. (1974), second-order slowly varying forces.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    query = mod(query, 2 * pi);
    extended = [theta;2 * pi];
    extended_values = [values;values(1)];
    value = interp1(extended, extended_values, query,'linear');
end
