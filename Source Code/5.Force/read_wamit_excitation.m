function reference = read_wamit_excitation(filename, rho, g, wave_amplitude, characteristic_length)
% READ_WAMIT_EXCITATION Read wamit excitation for the CRESTU hydrodynamic workflow.
%
% Syntax:
%   reference = read_wamit_excitation(filename, rho, g, wave_amplitude, characteristic_length)
%
% Description:
%   The routine converts first-order potential-flow or rigid-body data into generalized hydrodynamic coefficients, loads, restoring terms, or motions. Translational and rotational quantities retain the global 6-DOF ordering used by CRESTU.
%
% Inputs:
%   filename           - [character vector or string scalar] Input or output file path.
%   rho                - [scalar] Fluid density, [kg/m^3].
%   g                  - [scalar] Gravitational acceleration, [m/s^2].
%   wave_amplitude     - [scalar] Reference incident-wave amplitude, [m].
%   characteristic_length - [scalar] Characteristic length used for force and moment normalization, [m].
%
% Outputs:
%   reference          - [struct] Imported reference hydrodynamic data in documented SI normalization.
%
% Governing Equations / Theory:
%   Linear unsteady Bernoulli pressure, generalized surface integration, hydrostatics, radiation energy, or the frequency-domain rigid-body equation as applicable.
%
% References:
%   - Newman, J. N. (1977), Marine Hydrodynamics; Faltinsen, O. M. (1990), Sea Loads on Ships and Offshore Structures.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% --- 1. Validate Inputs and Initialize the Algorithm ---

    if nargin < 4 || isempty(wave_amplitude), wave_amplitude = 1; end
    if nargin < 5 || isempty(characteristic_length), characteristic_length = 1; end
    raw = readmatrix(filename, 'FileType', 'text');
    if size(raw, 2) < 7, error('CRESTU:WamitFormat', 'Expected seven columns in %s.', filename); end
    periods = unique(raw(:, 1), 'stable'); headings = unique(raw(:, 2), 'stable');
    nf = numel(periods); nh = numel(headings); ndof = max(raw(:, 3));
    value = complex(zeros(ndof, nh, nf)); coefficient = complex(zeros(ndof, nh, nf));
    for r = 1:size(raw, 1)
        k = find(periods == raw(r, 1), 1); h = find(headings == raw(r, 2), 1); mode = raw(r, 3);
        coefficient(mode, h, k) = complex(raw(r, 6), raw(r, 7));
        scale = rho * g * wave_amplitude * characteristic_length^(2 + (mode > 3));
        value(mode, h, k) = scale * coefficient(mode, h, k);
    end
    reference = struct('file', filename, 'periods', periods(:).', 'omegas', 2 * pi ./ periods(:).', ...
        'headings', headings(:).', 'coefficient', coefficient, 'force', value);
end
