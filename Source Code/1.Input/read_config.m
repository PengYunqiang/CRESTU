function cfg = read_config(config_file)
% READ_CONFIG Read config for the CRESTU hydrodynamic workflow.
%
% Syntax:
%   cfg = read_config(config_file)
%
% Description:
%   Parses and validates CRESTU configuration data.
%   Physical values and SI units are preserved.
%
% Inputs:
%   config_file        - [character vector or string scalar] CRESTU configuration-file path.
%
% Outputs:
%   cfg                - [struct] Parsed and validated CRESTU configuration in SI units.
%
% Governing Equations / Theory:
%   Deterministic configuration parsing and dimensional consistency checks.
%
% References:
%   - CRESTU project configuration specification.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    if nargin < 1 || isempty(config_file)
        config_file ='CRESTU.cfg';
    end
    config_file = absolute_path(config_file);
    fid = fopen(config_file,'rt');
    if fid < 0
        error('CRESTU:ConfigOpen','Cannot open configuration file: %s', config_file);
    end
    cleanup = onCleanup(@() fclose(fid));
    lines = cell(0, 1);
    while ~feof(fid)
        raw = fgetl(fid);
        if ~ischar(raw)
            continue;
        end
        raw = strtrim(normalize_text(raw));
        if isempty(raw) || startsWith(raw,'#') || startsWith(raw,'%')
            continue;
        end
        lines{end + 1, 1} = raw; %#ok<AGROW>
    end
    clear cleanup
    sections = split_sections(lines);
    required = {'PARA1','PARA2','PARA3','PARA4','PARA5','PARA6','PARA7','PARA8','PARA9','PARA10'};
    for k = 1:numel(required)
        if ~isfield(sections, required{k})
            error('CRESTU:MissingConfigSection','Missing %s in %s.', required{k}, config_file);
        end
    end

    cfg = struct('config_file', config_file,'config_dir', fileparts(config_file));
    cfg.case_name = strtrim(require_line(sections.PARA1, 1,'PARA1'));
    if isempty(cfg.case_name)
        error('CRESTU:InvalidCaseName','PARA1 case name is empty.');
    end

    numericValues = numeric_row(require_line(sections.PARA2, 1,'PARA2'),'PARA2');
    if numel(numericValues) ~= 4 && numel(numericValues) ~= 5
        error('CRESTU:ConfigArrayShape','PARA2 requires 4 flags, or 5 including IDRIFT.');
    end
    validate_flags(numericValues,'PARA2');
    idrift = 0;
    if numel(numericValues) == 5
        idrift = numericValues(5);
    end
    cfg.run = struct('ipoten', numericValues(1),'iforce', numericValues(2),'irad', numericValues(3),'idiff', numericValues(4), ...
'idrift', idrift);

    numericValues = numeric_row(require_line(sections.PARA3, 1,'PARA3'),'PARA3 frequency line');
    if numel(numericValues) ~= 4 && numel(numericValues) ~= 5
        error('CRESTU:ConfigArrayShape','PARA3 frequency line requires 4 or 5 values.');
    end
    cfg.freq.type = numericValues(1);
    cfg.freq.nfreq = numericValues(2);
    if cfg.freq.nfreq < 1 || fix(cfg.freq.nfreq) ~= cfg.freq.nfreq
        error('CRESTU:InvalidFrequencyCount','NFREQ must be a positive integer.');
    end
    base = linspace(numericValues(3), numericValues(4), cfg.freq.nfreq);
    expected_step = iff(cfg.freq.nfreq > 1, (numericValues(4) - numericValues(3)) / (cfg.freq.nfreq - 1), 0);
    cfg.freq.step = iff(numel(numericValues) == 5, numericValues(end), expected_step);
    if numel(numericValues) == 5 && cfg.freq.nfreq > 1 && abs(cfg.freq.step - expected_step) > 1e-10 * max(1, abs(expected_step))
        warning('CRESTU:FrequencyStepMismatch', ...
'FREQ_STEP (%g) conflicts with NFREQ/start/end (%g); NFREQ/start/end are authoritative.', ...
            cfg.freq.step, expected_step);
    end
    switch cfg.freq.type
        case 1
            cfg.freq.omegas = base;
        case 2
            if any(base <= 0)
                error('CRESTU:InvalidPeriods','Periods must be positive.');
            end
            cfg.freq.omegas = 2 * pi ./ base;
        case 3
            cfg.freq.omegas = 2 * pi .* base;

        otherwise
            error('CRESTU:InvalidFrequencyType', ...
'FREQ_TYPE must be 1, 2, or 3.');
    end
    if any(~isfinite(cfg.freq.omegas)|cfg.freq.omegas <= 0)
        error('CRESTU:InvalidFrequencies','All angular frequencies must be finite and positive.');
    end
    cfg.wave.ndir = scalar_integer(require_line(sections.PARA3, 2,'PARA3'),'NDIR');
    cfg.wave.headings = numeric_row(require_line(sections.PARA3, 3,'PARA3'),'wave headings');
    if numel(cfg.wave.headings) ~= cfg.wave.ndir
        error('CRESTU:HeadingCountMismatch','NDIR=%d but %d headings were supplied.', ...
            cfg.wave.ndir, numel(cfg.wave.headings));
    end

    numericValues = numeric_row(require_line(sections.PARA4, 1,'PARA4'),'PARA4');
    require_count(numericValues, 3,'PARA4');
    cfg.water_depth = numericValues(1);
    cfg.grav = numericValues(2);
    cfg.rho = numericValues(3);
    if cfg.water_depth < 0 || cfg.grav <= 0 || cfg.rho <= 0
        error('CRESTU:InvalidEnvironment','Depth must be nonnegative; gravity and density must be positive.');
    end

    cfg.n_bodies = scalar_integer(require_line(sections.PARA5, 1,'PARA5'),'NBODY');
    if cfg.n_bodies < 1
        error('CRESTU:InvalidBodyCount','NBODY must be at least one.');
    end
    if numel(sections.PARA5) < cfg.n_bodies + 1
        error('CRESTU:BodyConfigCount','PARA5 does not contain %d body records.', cfg.n_bodies);
    end
    cfg.bodies = repmat(struct('id', 0,'mesh_file','','pos', zeros(1, 3),'yaw', 0), cfg.n_bodies, 1);
    for bodyIndex = 1:cfg.n_bodies
        tokens = regexp(require_line(sections.PARA5, bodyIndex + 1,'PARA5'),'\s+','split');
        if numel(tokens) ~= 6
            error('CRESTU:BodyRecordShape','Body record %d must contain 6 tokens.', bodyIndex);
        end
        body_id = parse_scalar(tokens{1}, sprintf('body %d ID', bodyIndex));
        mesh_name = tokens{2};
        if strcmpi(mesh_name,'auto')
            if cfg.n_bodies == 1
                mesh_name = sprintf('%s_body.bmf', cfg.case_name);
            else
                mesh_name = sprintf('%s_body%d.bmf', cfg.case_name, body_id);
            end
        end
        cfg.bodies(bodyIndex).id = body_id;
        cfg.bodies(bodyIndex).mesh_file = resolve_path(cfg.config_dir, mesh_name);
        cfg.bodies(bodyIndex).pos = cellfun(@(s)parse_scalar(s,'body position'), tokens(3:5));
        cfg.bodies(bodyIndex).yaw = parse_scalar(tokens{6},'body yaw');
    end
    if numel(unique([cfg.bodies.id])) ~= cfg.n_bodies
        error('CRESTU:DuplicateBodyId','Body IDs must be unique.');
    end

    if numel(sections.PARA6) < 2 * cfg.n_bodies
        error('CRESTU:MassConfigCount','PARA6 requires two lines per body.');
    end
    cfg.mass_props = repmat(struct('mass', 0,'cg', zeros(1, 3),'inertia', zeros(3)), cfg.n_bodies, 1);
    for bodyIndex = 1:cfg.n_bodies
        m = numeric_row(sections.PARA6{2 * bodyIndex - 1}, sprintf('body %d mass/CG', bodyIndex));
        inertiaValues = numeric_row(sections.PARA6{2 * bodyIndex}, sprintf('body %d inertia', bodyIndex));
        require_count(m, 4,'mass and CG');
        require_count(inertiaValues, 6,'inertia tensor');
        cfg.mass_props(bodyIndex).mass = m(1);
        cfg.mass_props(bodyIndex).cg = m(2:4);
        cfg.mass_props(bodyIndex).inertia = [inertiaValues(1), -inertiaValues(4), -inertiaValues(5); -inertiaValues(4), inertiaValues(2), -inertiaValues(6); -inertiaValues(5), -inertiaValues(6), inertiaValues(3)];
        if cfg.mass_props(bodyIndex).mass <= 0 || any(eig(cfg.mass_props(bodyIndex).inertia) <= 0)
            error('CRESTU:InvalidMassProperties','Body %d mass/inertia must be positive definite.', bodyIndex);
        end
    end

    cfg.calc_modes = numeric_row(require_line(sections.PARA7, 1,'PARA7'),'PARA7 modes');
    if any(cfg.calc_modes < 1|cfg.calc_modes > 6|fix(cfg.calc_modes) ~= cfg.calc_modes)
        error('CRESTU:InvalidModes','Radiation modes must be integers in [1,6].');
    end
    numericValues = numeric_row(require_line(sections.PARA8, 1,'PARA8'),'PARA8');
    require_count(numericValues, 2,'PARA8');
    validate_flags(numericValues,'PARA8');
    cfg.isx = numericValues(1);
    cfg.isy = numericValues(2);
% is x>=0 for ISX and y>=0 for ISY.  Polar-vector and axial-vector
% parities follow reflection tensor rules.
    cfg.symmetry = struct();
    cfg.symmetry.multiplicity = 2^(cfg.isx + cfg.isy);
    cfg.symmetry.area_scale = cfg.symmetry.multiplicity;
    cfg.symmetry.x_bounds = iff(cfg.isx == 1, [0, Inf], [-Inf, Inf]);
    cfg.symmetry.y_bounds = iff(cfg.isy == 1, [0, Inf], [-Inf, Inf]);
    cfg.symmetry.mode_parity_x = [-1, 1, 1, 1, -1, -1]; % surge,sway,heave,roll,pitch,yaw
    cfg.symmetry.mode_parity_y = [1, -1, 1, -1, 1, -1];

    numericValues = numeric_row(require_line(sections.PARA9, 1,'PARA9'),'PARA9');
    require_count(numericValues, 7,'PARA9');
    cfg.z_tol = numericValues(1); cfg.fs = struct('nr_near', numericValues(2),'nr_sponge', numericValues(3),'r_inner', numericValues(4), ...
'r_outer', numericValues(5),'sponge_ratio', numericValues(6),'nz_farfield', numericValues(7));
    cfg.fs.mu0 = parse_scalar(require_line(sections.PARA10, 1,'PARA10'),'PARA10');
    if cfg.z_tol <= 0 || any([cfg.fs.nr_near, cfg.fs.nr_sponge, cfg.fs.nz_farfield] < 1) || ...
            cfg.fs.r_inner <= 0 || cfg.fs.r_outer <= cfg.fs.r_inner || ...
            cfg.fs.sponge_ratio <= 0 || cfg.fs.mu0 < 0
        error('CRESTU:InvalidDomainParameters','PARA9/PARA10 parameters are invalid.');
    end

    cfg.files.fs = fullfile(cfg.config_dir, sprintf('%s_fs.bmf', cfg.case_name));
    cfg.files.seabed = fullfile(cfg.config_dir, sprintf('%s_seabed.bmf', cfg.case_name));
    cfg.files.farfield = fullfile(cfg.config_dir, sprintf('%s_farfield.bmf', cfg.case_name));
    cfg.files.potential_cache = fullfile(cfg.config_dir, sprintf('%s_PotCache.mat', cfg.case_name));
    cfg.files.results = fullfile(cfg.config_dir, sprintf('%s_Results.mat', cfg.case_name));
    cfg.outer_reference_mesh_file = parse_outer_reference_mesh(sections, cfg.config_dir);
    fprintf('[OK] Loaded %s\n', config_file);
    fprintf('    Case=%s, IPOTEN=%d, IFORCE=%d, IRAD=%d, IDIFF=%d, IDRIFT=%d\n', cfg.case_name, ...
        cfg.run.ipoten, cfg.run.iforce, cfg.run.irad, cfg.run.idiff, cfg.run.idrift);
    fprintf('    Bodies=%d, frequencies=%d, headings=%d, depth=%.3g m\n', ...
        cfg.n_bodies, cfg.freq.nfreq, cfg.wave.ndir, cfg.water_depth);
end

function meshFilename = parse_outer_reference_mesh(sections, configDirectory)
% PARSE_OUTER_REFERENCE_MESH Read an optional fixed waterline reference mesh.
%
% Syntax:
%   meshFilename = parse_outer_reference_mesh(sections, configDirectory)
%
% Inputs:
%   sections        - Parsed configuration sections [-].
%   configDirectory - Configuration-file directory [-].
%
% Outputs:
%   meshFilename    - Absolute reference BMF path, or empty text [-].

    %% 阶段 1: 处理未配置的可选段

    meshFilename = '';

    if ~isfield(sections, 'PARA11')
        return
    end

    assert(numel(sections.PARA11) == 1, ...
        'CRESTU:OuterReferenceShape', ...
        'PARA11 requires exactly one OUTER_WATERLINE_MESH record.');

    %% 阶段 2: 解析并验证外域水线参考网格

    tokens = regexp(strtrim(sections.PARA11{1}), '\s+', 'split');
    assert(numel(tokens) == 2 && ...
        strcmpi(tokens{1}, 'OUTER_WATERLINE_MESH'), ...
        'CRESTU:OuterReferenceSyntax', ...
        'PARA11 syntax is: OUTER_WATERLINE_MESH path-or-NONE.');

    if strcmpi(tokens{2}, 'NONE')
        return
    end

    meshFilename = resolve_path(configDirectory, tokens{2});
    assert(isfile(meshFilename), 'CRESTU:OuterReferenceMissing', ...
        'Outer-domain waterline reference mesh was not found: %s', ...
        meshFilename);
end

function sections = split_sections(lines)
% SPLIT_SECTIONS Partition configuration text into named sections.
%
% Syntax:
%   sections = split_sections(lines)
%
% Description:
%   Parses and validates CRESTU configuration data.
%   Physical values and SI units are preserved.
%
% Inputs:
%   lines              - [cell array] Configuration-file lines.
%
% Outputs:
%   sections           - [struct] Configuration lines partitioned by section.
%
% Governing Equations / Theory:
%   Deterministic configuration parsing and dimensional consistency checks.
%
% References:
%   - CRESTU project configuration specification.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    sections = struct();
    current ='';
    for k = 1:numel(lines)
        token = regexp(lines{k},'^\s*(PARA\d+)\s*:','tokens','once','ignorecase');
        if ~isempty(token)
            current = upper(token{1});
            sections.(current) = cell(0, 1);
        elseif isempty(current)
            error('CRESTU:ConfigPreamble','Data appears before a PARA section: %s', lines{k});
        else
            values = sections.(current);
            values{end + 1, 1} = lines{k}; %#ok<AGROW> section lengths are unknown while parsing
            sections.(current) = values;
        end
    end
end
function line = require_line(lines, idx, section)
% REQUIRE_LINE Return a required configuration line or raise a validation error.
%
% Syntax:
%   line = require_line(lines, idx, section)
%
% Description:
%   Parses and validates CRESTU configuration data.
%   Physical values and SI units are preserved.
%
% Inputs:
%   lines              - [cell array] Configuration-file lines.
%   idx                - [scalar] One-based line index, dimensionless.
%   section            - [character vector] Configuration section name.
%
% Outputs:
%   line               - [character vector] Validated configuration line.
%
% Governing Equations / Theory:
%   Deterministic configuration parsing and dimensional consistency checks.
%
% References:
%   - CRESTU project configuration specification.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    if numel(lines) < idx
        error('CRESTU:MissingConfigData','%s is missing line %d.', section, idx);
    end
    line = lines{idx};
end
function values = numeric_row(line, label)
% NUMERIC_ROW Parse one configuration row as finite numeric values.
%
% Syntax:
%   values = numeric_row(line, label)
%
% Description:
%   Parses and validates CRESTU configuration data.
%   Physical values and SI units are preserved.
%
% Inputs:
%   line               - [character vector] Configuration line to parse.
%   label              - [character vector or string scalar] Diagnostic field label.
%
% Outputs:
%   values             - [M x K] Interpolated samples with the same physical units as the input field.
%
% Governing Equations / Theory:
%   Deterministic configuration parsing and dimensional consistency checks.
%
% References:
%   - CRESTU project configuration specification.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    tokens = regexp(strtrim(normalize_text(line)),'[,\s]+','split');
    values = zeros(1, numel(tokens));
    for k = 1:numel(tokens)
        values(k) = parse_scalar(tokens{k}, label);
    end
end
function value = parse_scalar(token, label)
% PARSE_SCALAR Parse and validate one scalar configuration token.
%
% Syntax:
%   value = parse_scalar(token, label)
%
% Description:
%   Parses and validates CRESTU configuration data.
%   Physical values and SI units are preserved.
%
% Inputs:
%   token              - [character vector or string scalar] Scalar configuration token.
%   label              - [character vector or string scalar] Diagnostic field label.
%
% Outputs:
%   value              - [scalar] Interpolated or parsed value with units inherited from the input contract.
%
% Governing Equations / Theory:
%   Deterministic configuration parsing and dimensional consistency checks.
%
% References:
%   - CRESTU project configuration specification.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    token = normalize_text(token);
    value = str2double(token);
    if ~isscalar(value) || ~isfinite(value)
        error('CRESTU:InvalidNumber','Invalid token "%s" in %s.', token, label);
    end
end
function value = scalar_integer(line, label)
% SCALAR_INTEGER Parse and validate one integer-valued configuration field.
%
% Syntax:
%   value = scalar_integer(line, label)
%
% Description:
%   Parses and validates CRESTU configuration data.
%   Physical values and SI units are preserved.
%
% Inputs:
%   line               - [character vector] Configuration line to parse.
%   label              - [character vector or string scalar] Diagnostic field label.
%
% Outputs:
%   value              - [scalar] Interpolated or parsed value with units inherited from the input contract.
%
% Governing Equations / Theory:
%   Deterministic configuration parsing and dimensional consistency checks.
%
% References:
%   - CRESTU project configuration specification.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    values = numeric_row(line, label);
    require_count(values, 1, label);
    value = values(1);
    if value < 0 || fix(value) ~= value
        error('CRESTU:InvalidInteger','%s must be a nonnegative integer.', label);
    end
end
function require_count(values, count, label)
% REQUIRE_COUNT Validate the number of values in a configuration row.
%
% Syntax:
%   require_count(values, count, label)
%
% Description:
%   Parses and validates CRESTU configuration data.
%   Physical values and SI units are preserved.
%
% Inputs:
%   values             - [numeric array] Samples to be transformed or interpolated; units are preserved.
%   count              - [scalar] Required value count, dimensionless.
%   label              - [character vector or string scalar] Diagnostic field label.
%
% Outputs:
%   None.
%
% Governing Equations / Theory:
%   Deterministic configuration parsing and dimensional consistency checks.
%
% References:
%   - CRESTU project configuration specification.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    if numel(values) ~= count
        error('CRESTU:ConfigArrayShape','%s requires %d values; got %d.', label, count, numel(values));
    end
end
function validate_flags(values, label)
% VALIDATE_FLAGS Validate binary configuration flags.
%
% Syntax:
%   validate_flags(values, label)
%
% Description:
%   Parses and validates CRESTU configuration data.
%   Physical values and SI units are preserved.
%
% Inputs:
%   values             - [numeric array] Samples to be transformed or interpolated; units are preserved.
%   label              - [character vector or string scalar] Diagnostic field label.
%
% Outputs:
%   None.
%
% Governing Equations / Theory:
%   Deterministic configuration parsing and dimensional consistency checks.
%
% References:
%   - CRESTU project configuration specification.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    if any(values ~= 0&values ~= 1)
        error('CRESTU:InvalidFlag','%s flags must be 0 or 1.', label);
    end
end
function out = normalize_text(in)
% NORMALIZE_TEXT Normalize a configuration value to a trimmed character vector.
%
% Syntax:
%   out = normalize_text(in)
%
% Description:
%   Parses and validates CRESTU configuration data.
%   Physical values and SI units are preserved.
%
% Inputs:
%   in                 - [character vector, string scalar, or numeric value] Value to normalize.
%
% Outputs:
%   out                - [value] Normalized or selected result with type and physical units preserved from the input.
%
% Governing Equations / Theory:
%   Deterministic configuration parsing and dimensional consistency checks.
%
% References:
%   - CRESTU project configuration specification.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    out = strrep(in, char(8722),'-');
    out = strrep(out, char(8208),'-');
    out = strrep(out, char(8209),'-');
    out = strrep(out, char(8211),'-');
    out = strrep(out, char(160),' ');
end
function out = absolute_path(in)
% ABSOLUTE_PATH Convert an input path to an absolute normalized path.
%
% Syntax:
%   out = absolute_path(in)
%
% Description:
%   Parses and validates CRESTU configuration data.
%   Physical values and SI units are preserved.
%
% Inputs:
%   in                 - [character vector, string scalar, or numeric value] Value to normalize.
%
% Outputs:
%   out                - [value] Normalized or selected result with type and physical units preserved from the input.
%
% Governing Equations / Theory:
%   Deterministic configuration parsing and dimensional consistency checks.
%
% References:
%   - CRESTU project configuration specification.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    if is_absolute(in)
        out = in;
    else
        out = fullfile(pwd, in);
    end
end
function out = resolve_path(base_dir, in)
% RESOLVE_PATH Resolve a project-relative path against a base directory.
%
% Syntax:
%   out = resolve_path(base_dir, in)
%
% Description:
%   Parses and validates CRESTU configuration data.
%   Physical values and SI units are preserved.
%
% Inputs:
%   base_dir           - [character vector or string scalar] Absolute base-directory path.
%   in                 - [character vector, string scalar, or numeric value] Value to normalize.
%
% Outputs:
%   out                - [value] Normalized or selected result with type and physical units preserved from the input.
%
% Governing Equations / Theory:
%   Deterministic configuration parsing and dimensional consistency checks.
%
% References:
%   - CRESTU project configuration specification.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    if is_absolute(in)
        out = in;
    else
        out = fullfile(base_dir, in);
    end
end
function tf = is_absolute(in)
% IS_ABSOLUTE Test whether a path is absolute on the active platform.
%
% Syntax:
%   tf = is_absolute(in)
%
% Description:
%   Parses and validates CRESTU configuration data.
%   Physical values and SI units are preserved.
%
% Inputs:
%   in                 - [character vector, string scalar, or numeric value] Value to normalize.
%
% Outputs:
%   tf                 - [logical scalar] True when the tested condition is satisfied.
%
% Governing Equations / Theory:
%   Deterministic configuration parsing and dimensional consistency checks.
%
% References:
%   - CRESTU project configuration specification.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    tf = ~isempty(regexp(in,'^[A-Za-z]:[\\/]','once')) || startsWith(in, filesep);
end
function out = iff(condition, a, bodyIndex)
% IFF Select one of two values from a logical condition.
%
% Syntax:
%   out = iff(condition, a, bodyIndex)
%
% Description:
%   Parses and validates CRESTU configuration data.
%   Physical values and SI units are preserved.
%
% Inputs:
%   condition          - [logical scalar] Selection condition.
%   a                  - First vector or scalar operand; dimensions and SI units follow the stated algorithm.
%   bodyIndex                  - Second vector or scalar operand; dimensions and SI units follow the stated algorithm.
%
% Outputs:
%   out                - [value] Normalized or selected result with type and physical units preserved from the input.
%
% Governing Equations / Theory:
%   Deterministic configuration parsing and dimensional consistency checks.
%
% References:
%   - CRESTU project configuration specification.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    if condition
        out = a;
    else
        out = bodyIndex;
    end
end
