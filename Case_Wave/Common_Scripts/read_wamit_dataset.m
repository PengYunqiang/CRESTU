function dataset = read_wamit_dataset(folder, characteristic_length)
% READ_WAMIT_DATASET Read a complete first- and second-order WAMIT dataset.
%
% Syntax:
%   dataset = read_wamit_dataset(folder, characteristic_length)
%
% Inputs:
%   folder                : [char|string] Read-only WAMIT result directory.
%   characteristic_length : [scalar] CRESTU reporting length, in m.
%
% Outputs:
%   dataset               : [struct] Metadata and parsed .1/.2/.3/.4/.8/.9 records.
%
% Mathematical Reference:
%   WAMIT User Manual output-file conventions and dimensional scalings.
    if nargin < 2 || isempty(characteristic_length)
        characteristic_length = 1;
    end
    folder = char(folder);
    if ~isfolder(folder)
        error('CRESTU:WamitFolder', 'WAMIT directory does not exist: %s', folder);
    end

% ==========================================
% Read metadata from the text report
% ==========================================
    metadata = struct('gravity', 9.80665, 'density', 1025, 'water_depth', 0, ...
        'length_scale', 1, 'body_count', NaN, 'body_positions', zeros(0, 3));
    out_file = find_extension_file(folder, 'out', false);
    if ~isempty(out_file)
        report_text = fileread(out_file);
        metadata.gravity = extract_scalar(report_text, 'Gravity:\s*([+\-0-9.Ee]+)', metadata.gravity);
        metadata.density = extract_scalar(report_text, 'Water density:\s*([+\-0-9.Ee]+)', metadata.density);
        metadata.water_depth = extract_scalar(report_text, 'Water depth:\s*([+\-0-9.Ee]+)', metadata.water_depth);
        metadata.length_scale = extract_scalar(report_text, 'Length scale:\s*([+\-0-9.Ee]+)', ...
            metadata.length_scale);
        position_tokens = regexp(report_text, ...
            'XBODY\s*=\s*([+\-0-9.Ee]+)\s+YBODY\s*=\s*([+\-0-9.Ee]+)\s+ZBODY\s*=\s*([+\-0-9.Ee]+)', ...
            'tokens');
        if ~isempty(position_tokens)
            metadata.body_positions = cellfun(@str2double, vertcat(position_tokens{:}));
            metadata.body_count = size(metadata.body_positions, 1);
        end
    end

% ==========================================
% Parse all available hydrodynamic outputs
% ==========================================
    dataset = struct('folder', folder, 'metadata', metadata, 'radiation', [], ...
        'excitation', [], 'scattering', [], 'rao', [], 'drift_momentum', [], ...
        'drift_pressure', [], 'files', struct());

    dataset.files.one = find_extension_file(folder, '1', true);
    dataset.radiation = read_wamit_first_order(dataset.files.one, metadata.density);

    dataset.files.two = find_extension_file(folder, '2', false);
    if ~isempty(dataset.files.two)
        dataset.excitation = read_wamit_excitation(dataset.files.two, metadata.density, ...
            metadata.gravity, 1, metadata.length_scale);
    end
    dataset.files.three = find_extension_file(folder, '3', false);
    if ~isempty(dataset.files.three)
        dataset.scattering = read_wamit_excitation(dataset.files.three, metadata.density, ...
            metadata.gravity, 1, metadata.length_scale);
    end
    dataset.files.four = find_extension_file(folder, '4', false);
    if ~isempty(dataset.files.four)
        dataset.rao = read_wamit_rao(dataset.files.four);
    end
    dataset.files.eight = find_extension_file(folder, '8', false);
    if ~isempty(dataset.files.eight)
        dataset.drift_momentum = read_wamit_mean_drift(dataset.files.eight, metadata.density, ...
            metadata.gravity, 1, metadata.length_scale, characteristic_length);
    end
    dataset.files.nine = find_extension_file(folder, '9', false);
    if ~isempty(dataset.files.nine)
        dataset.drift_pressure = read_wamit_mean_drift(dataset.files.nine, metadata.density, ...
            metadata.gravity, 1, metadata.length_scale, characteristic_length);
    end
end

function filename = find_extension_file(folder, extension, required)
% FIND_EXTENSION_FILE Locate the unique WAMIT file for an extension.
%
% Syntax:
%   filename = find_extension_file(folder, extension, required)
%
% Inputs:
%   folder    : [char] WAMIT result directory.
%   extension : [char] Extension without a leading dot.
%   required  : [logical scalar] Throw if no matching file exists.
%
% Outputs:
%   filename  : [char] Full path or an empty character vector.
%
% Mathematical Reference:
%   File-system utility; no mathematical model is used.
    candidates = dir(fullfile(folder, ['*.' extension]));
    candidates = candidates(~[candidates.isdir]);
    if isempty(candidates)
        if required
            error('CRESTU:WamitFile', 'No .%s file exists in %s.', extension, folder);
        end
        filename = '';
        return
    end
    [~, selected] = min([candidates.bytes]);
    filename = fullfile(candidates(selected).folder, candidates(selected).name);
end

function value = extract_scalar(text, expression, default_value)
% EXTRACT_SCALAR Extract one floating-point value with a regular expression.
%
% Syntax:
%   value = extract_scalar(text, expression, default_value)
%
% Inputs:
%   text          : [char] Source text.
%   expression    : [char] Regular expression with one capture group.
%   default_value : [scalar] Value returned when the expression is absent.
%
% Outputs:
%   value         : [scalar] Parsed finite value or the supplied default.
%
% Mathematical Reference:
%   Text parsing utility; no mathematical model is used.
    token = regexp(text, expression, 'tokens', 'once');
    if isempty(token)
        value = default_value;
    else
        value = str2double(token{1});
    end
end
