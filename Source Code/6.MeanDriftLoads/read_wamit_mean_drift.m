function reference = read_wamit_mean_drift(filename, rho, g, wave_amplitude, wamit_length, output_length)
% READ_WAMIT_MEAN_DRIFT Read wamit mean drift for the CRESTU hydrodynamic workflow.
%
% Syntax:
%   reference = read_wamit_mean_drift(filename, rho, g, wave_amplitude, wamit_length, output_length)
%
% Description:
%   Computes or processes second-order mean wave-drift quantities.
%   Complex averages follow exp(i*omega*t) and the 6-DOF order.
%
% Inputs:
%   filename           - [character vector or string scalar] Input or output file path.
%   rho                - [scalar] Fluid density, [kg/m^3].
%   g                  - [scalar] Gravitational acceleration, [m/s^2].
%   wave_amplitude     - [scalar] Reference incident-wave amplitude, [m].
%   wamit_length       - [scalar] Characteristic length used by WAMIT normalization, [m].
%   output_length      - [scalar] Characteristic length required for output normalization, [m].
%
% Outputs:
%   reference          - [struct] Imported reference hydrodynamic data in documented SI normalization.
%
% Governing Equations / Theory:
%   Pinkster near-field direct pressure integration, Maruo-Newman far-field momentum balance, or supporting surface reconstruction as applicable.
%
% References:
%   - Pinkster, J. A. (1980), Low Frequency Second Order Wave Exciting Forces on Floating Structures; Newman, J. N. (1974), second-order slowly varying forces.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    if nargin < 4 || isempty(wave_amplitude)
        wave_amplitude = 1;
    end
    if nargin < 5 || isempty(wamit_length)
        wamit_length = 1;
    end
    if nargin < 6 || isempty(output_length)
        output_length = wamit_length;
    end

    raw = readmatrix(filename,'FileType','text');
    if size(raw, 2) < 8
        error('CRESTU:WamitFormat','Expected eight columns in %s.', filename);
    end
    raw = raw(raw(:, 2) == raw(:, 3) & raw(:, 4) > 0 & raw(:, 4) <= 6, :);
    periods = unique(raw(:, 1),'stable');
    headings = unique(raw(:, 2),'stable');
    frequency_count = numel(periods);
    heading_count = numel(headings);
    coefficient = zeros(6, heading_count, frequency_count);
    load_value = zeros(6, heading_count, frequency_count);
    cd = zeros(6, heading_count, frequency_count);

    for row_index = 1:size(raw, 1)
        frequency_index = find(periods == raw(row_index, 1), 1);
        heading_index = find(headings == raw(row_index, 2), 1);
        mode = raw(row_index, 4);
        coefficient(mode, heading_index, frequency_index) = raw(row_index, 7);
        wamit_scale = rho * g * wave_amplitude^2 * wamit_length^(1 + (mode > 3));
        output_scale = 0.5 * rho * g * wave_amplitude^2 * output_length^(1 + (mode > 3));
        load_value(mode, heading_index, frequency_index) = wamit_scale * raw(row_index, 7);
        cd(mode, heading_index, frequency_index) = ...
            load_value(mode, heading_index, frequency_index) / output_scale;
    end

    [~, ~, extension] = fileparts(filename);
    method ='direct_pressure';
    if strcmpi(extension,'.8')
        method ='momentum';
    end
    reference = struct( ...
'file', filename, ...
'method', method, ...
'periods', periods(:).', ...
'omegas', 2 * pi ./ periods(:).', ...
'headings', headings(:).', ...
'wamit_coefficient', coefficient, ...
'Cd', cd, ...
'load', load_value, ...
'wamit_length', wamit_length, ...
'output_length', output_length);
end
