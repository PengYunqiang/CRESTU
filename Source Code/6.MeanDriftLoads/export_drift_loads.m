function output = export_drift_loads(filename, omegas, headings, near_loads, far_loads, rho, g, wave_amplitude, L)
% EXPORT_DRIFT_LOADS Export drift loads for the CRESTU hydrodynamic workflow.
%
% Syntax:
%   output = export_drift_loads(filename, omegas, headings, near_loads, far_loads, rho, g, wave_amplitude, L)
%
% Description:
%   The routine evaluates, reconstructs, imports, or exports quantities required by second-order mean wave-drift analysis. Complex products are time averaged consistently with the exp(i*omega*t) convention and generalized loads use the project 6-DOF ordering.
%
% Inputs:
%   filename           - [character vector or string scalar] Input or output file path.
%   omegas             - [1 x Nf] Angular-frequency samples, [rad/s].
%   headings           - [1 x Nh] Wave headings, [deg].
%   near_loads         - [6Nb x Nh x Nf] Near-field mean generalized loads, [N] and [N m].
%   far_loads          - [6 x Nh x Nf] Global far-field mean loads, [N] and [N m].
%   rho                - [scalar] Fluid density, [kg/m^3].
%   g                  - [scalar] Gravitational acceleration, [m/s^2].
%   wave_amplitude     - [scalar] Reference incident-wave amplitude, [m].
%   L                  - [scalar] Characteristic body length, [m].
%
% Outputs:
%   output             - [struct] Exported load data and normalization metadata.
%
% Governing Equations / Theory:
%   Pinkster near-field direct pressure integration, Maruo-Newman far-field momentum balance, or supporting surface reconstruction as applicable.
%
% References:
%   - Pinkster, J. A. (1980), Low Frequency Second Order Wave Exciting Forces on Floating Structures; Newman, J. N. (1974), second-order slowly varying forces.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% --- 1. Validate Inputs and Initialize the Algorithm ---

    force_scale = 0.5 * rho * g * wave_amplitude^2 * L;
    moment_scale = force_scale * L;
    nf = numel(omegas); nh = numel(headings); ndof = size(near_loads, 1);
    if mod(ndof, 6) ~= 0 || size(far_loads, 1) ~= 6
        error('CRESTU:DriftExportShape', 'Near loads require 6N rows and global far-field loads require 6 rows.');
    end
    rows = nf * nh * (ndof + 6); frequency = zeros(rows, 1); heading = zeros(rows, 1);
    method = strings(rows, 1); body = zeros(rows, 1); dof = zeros(rows, 1); load_value = zeros(rows, 1); cursor = 0;
    for k = 1:nf
        for h = 1:nh
            idx = cursor + (1:ndof); frequency(idx) = omegas(k); heading(idx) = headings(h);
            method(idx) = "near_field"; body(idx) = repelem((1:ndof / 6).', 6); dof(idx) = repmat((1:6).', ndof / 6, 1);
            load_value(idx) = near_loads(:, h, k); cursor = cursor + ndof;
            idx = cursor + (1:6); frequency(idx) = omegas(k); heading(idx) = headings(h);
            method(idx) = "far_field_global"; body(idx) = 0; dof(idx) = (1:6).';
            load_value(idx) = far_loads(:, h, k); cursor = cursor + 6;
        end
    end
    scale = force_scale * ones(rows, 1); scale(dof > 3) = moment_scale;
    output = table(frequency, heading, method, body, dof, load_value, load_value ./ scale, ...
        'VariableNames', {'omega_rad_s', 'heading_deg', 'method', 'body', 'dof', 'load', 'Cd'});
    [folder, ~, ext] = fileparts(filename); if ~isempty(folder) && ~exist(folder, 'dir'), mkdir(folder); end
    switch lower(ext)
        case '.csv', writetable(output, filename);
        case '.mat', save(filename, 'output', 'force_scale', 'moment_scale');
        otherwise, error('CRESTU:DriftExportFormat', 'Use .csv or .mat for drift-load export.');
    end
    fprintf('>>> Mean drift loads exported: %s\n', filename);
end
