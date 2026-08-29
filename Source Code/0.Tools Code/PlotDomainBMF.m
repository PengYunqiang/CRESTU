clear;
clc;
close all;

%% =========================================================================
%% Stage 1: Select a File and Identify the Case Prefix
%% =========================================================================
[filename, filepath] = uigetfile('*.bmf;*.BMF','Select any BMF mesh file from the case');
if isequal(filename, 0)
    disp('The user canceled file selection.');
    return;
end

[~, raw_name, ~] = fileparts(filename);

% Strip component suffixes and identify the case name automatically
case_name = regexprep(raw_name,'(_fs|_seabed|_farfield|_body\d*|_body)$','','ignorecase');
fprintf('[INFO] Identified case prefix (Case Name): %s\n', case_name);

%% =========================================================================
%% Stage 2: Read Every BMF Boundary Component
%% =========================================================================
% 2.1 body-surface mesh (support either case.bmf or multibody case_body*.bmf files)
body_files = dir(fullfile(filepath, [case_name,'_body*.bmf']));
if isempty(body_files)
    single_body = fullfile(filepath, [case_name,'.bmf']);
    if exist(single_body,'file')
        body_files = dir(single_body);
    end
end

body_list = {};
if ~isempty(body_files)
    for bodyFileIndex = 1:length(body_files)
        b_path = fullfile(filepath, body_files(bodyFileIndex).name);
        mesh_b = read_bmf(b_path);

% Evaluate normals and centers when read_bmf does not supply them
        if ~isfield(mesh_b,'centers') || ~isfield(mesh_b,'normals')
            [mesh_b.centers, mesh_b.normals] = compute_normals_and_centers(mesh_b.vertices);
        end

        body_list{end + 1} = mesh_b; %#ok<SAGROW>
        fprintf('  [OK] Loaded body surface: %s (NPAN=%d)\n', body_files(bodyFileIndex).name, mesh_b.n_panels);
    end
else
    warning('Could not find a body-surface mesh (%s.bmf or %s_body*.bmf).', case_name, case_name);
end

% 2.2 free surface (Free Surface)
file_fs = fullfile(filepath, sprintf('%s_fs.bmf', case_name));
mesh_fs = [];
if exist(file_fs,'file')
    mesh_fs = read_bmf(file_fs);
    if ~isfield(mesh_fs,'mu_damping')
        mesh_fs.mu_damping = zeros(1, mesh_fs.n_panels); % set the default damping value to zero
    end
    fprintf('  [OK] Loadedfree surface: %s_fs.bmf (NPAN=%d)\n', case_name, mesh_fs.n_panels);
else
    warning('Free-surface mesh is missing: %s', file_fs);
end

% 2.3 seabed mesh (Seabed)
file_seabed = fullfile(filepath, sprintf('%s_seabed.bmf', case_name));
mesh_seabed = [];
if exist(file_seabed,'file')
    mesh_seabed = read_bmf(file_seabed);
    fprintf('  [OK] Loadedseabed mesh: %s_seabed.bmf (NPAN=%d)\n', case_name, mesh_seabed.n_panels);
else
    warning('seabed mesh is missing: %s', file_seabed);
end

% 2.4 far-field cylinder (Farfield)
file_farfield = fullfile(filepath, sprintf('%s_farfield.bmf', case_name));
mesh_farfield = [];
if exist(file_farfield,'file')
    mesh_farfield = read_bmf(file_farfield);
    fprintf('  [OK] Loadedfar-field mesh: %s_farfield.bmf (NPAN=%d)\n', case_name, mesh_farfield.n_panels);
else
    warning('far-field mesh is missing: %s', file_farfield);
end

%% =========================================================================
%% Stage 3: Assemble the Domain and Infer Missing Dimensions
%% =========================================================================
% Estimate the outer radius r_outer
r_outer = 50.0; % default fallback value
if ~isempty(mesh_fs)
    r_outer = max(sqrt(mesh_fs.vertices(:, :, 1) .^ 2 + mesh_fs.vertices(:, :, 2) .^ 2), [],'all');
elseif ~isempty(mesh_farfield)
    r_outer = max(sqrt(mesh_farfield.vertices(:, :, 1) .^ 2 + mesh_farfield.vertices(:, :, 2) .^ 2), [],'all');
end

% Estimate the water depth water_depth
water_depth = 0;
if ~isempty(mesh_seabed)
    water_depth = abs(mean(mesh_seabed.vertices(:, :, 3),'all'));
end

% Construct a lightweight domain
cfg = struct();
cfg.case_name = case_name;
cfg.n_bodies = length(body_list);
cfg.water_depth = water_depth;
cfg.fs.r_outer = r_outer;

domain = struct();
domain.cfg = cfg;
domain.body_list = body_list;
domain.fs = mesh_fs;
domain.seabed = mesh_seabed;
domain.farfield = mesh_farfield;
domain.waterline = {}; % Leave empty when no explicit waterline is available

%% =========================================================================
%% Stage 4: Render the Complete Domain
%% =========================================================================
fig = plot_bmf_domain(domain,'full');