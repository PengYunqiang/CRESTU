clear; clc; close all;

%% =========================================================================
%% 1. 选择文件并提取 Case 前缀
%% =========================================================================
[filename, filepath] = uigetfile('*.bmf;*.BMF', '选择该算例下的任意一个 BMF 网格文件');
if isequal(filename, 0)
    disp('用户取消了文件选择。');
    return;
end

[~, raw_name, ~] = fileparts(filename);

% 剥离后缀 (_fs, _seabed, _farfield, _body, _body1 等)，自动提取算例名
case_name = regexprep(raw_name, '(_fs|_seabed|_farfield|_body\d*|_body)$', '', 'ignorecase');
fprintf('>>> 识别到算例前缀 (Case Name): %s\n', case_name);

%% =========================================================================
%% 2. 自动读取各个区域的 BMF 网格
%% =========================================================================
% 2.1 物面网格 (支持单体 case.bmf 或多体 case_body*.bmf)
body_files = dir(fullfile(filepath, [case_name, '_body*.bmf']));
if isempty(body_files)
    single_body = fullfile(filepath, [case_name, '.bmf']);
    if exist(single_body, 'file')
        body_files = dir(single_body);
    end
end

body_list = {};
if ~isempty(body_files)
    for b = 1:length(body_files)
        b_path = fullfile(filepath, body_files(b).name);
        mesh_b = read_bmf(b_path);
        
        % 补充法向与中心点（若 read_bmf 未直接计算）
        if ~isfield(mesh_b, 'centers') || ~isfield(mesh_b, 'normals')
            [mesh_b.centers, mesh_b.normals] = compute_normals_and_centers(mesh_b.vertices);
        end
        
        body_list{end+1} = mesh_b; %#ok<SAGROW>
        fprintf('  [√] 已加载物面: %s (NPAN=%d)\n', body_files(b).name, mesh_b.n_panels);
    end
else
    warning('未找到物面网格文件 (%s.bmf 或 %s_body*.bmf)。', case_name, case_name);
end

% 2.2 自由水面 (Free Surface)
file_fs = fullfile(filepath, sprintf('%s_fs.bmf', case_name));
mesh_fs = [];
if exist(file_fs, 'file')
    mesh_fs = read_bmf(file_fs);
    if ~isfield(mesh_fs, 'mu_damping')
        mesh_fs.mu_damping = zeros(1, mesh_fs.n_panels); % 缺省阻尼值置 0
    end
    fprintf('  [√] 已加载自由水面: %s_fs.bmf (NPAN=%d)\n', case_name, mesh_fs.n_panels);
else
    warning('自由水面网格缺失: %s', file_fs);
end

% 2.3 海底网格 (Seabed)
file_seabed = fullfile(filepath, sprintf('%s_seabed.bmf', case_name));
mesh_seabed = [];
if exist(file_seabed, 'file')
    mesh_seabed = read_bmf(file_seabed);
    fprintf('  [√] 已加载海底网格: %s_seabed.bmf (NPAN=%d)\n', case_name, mesh_seabed.n_panels);
else
    warning('海底网格缺失: %s', file_seabed);
end

% 2.4 远场柱面 (Farfield)
file_farfield = fullfile(filepath, sprintf('%s_farfield.bmf', case_name));
mesh_farfield = [];
if exist(file_farfield, 'file')
    mesh_farfield = read_bmf(file_farfield);
    fprintf('  [√] 已加载远场网格: %s_farfield.bmf (NPAN=%d)\n', case_name, mesh_farfield.n_panels);
else
    warning('远场网格缺失: %s', file_farfield);
end

%% =========================================================================
%% 3. 组装 domain 结构体 (自动推断缺失的尺寸参数)
%% =========================================================================
% 自动估算外半径 r_outer
r_outer = 50.0; % 缺省回退值
if ~isempty(mesh_fs)
    r_outer = max(sqrt(mesh_fs.vertices(:,:,1).^2 + mesh_fs.vertices(:,:,2).^2), [], 'all');
elseif ~isempty(mesh_farfield)
    r_outer = max(sqrt(mesh_farfield.vertices(:,:,1).^2 + mesh_farfield.vertices(:,:,2).^2), [], 'all');
end

% 自动估算水深 water_depth
water_depth = 0;
if ~isempty(mesh_seabed)
    water_depth = abs(mean(mesh_seabed.vertices(:,:,3), 'all'));
end

% 构建轻量化 domain
cfg = struct();
cfg.case_name   = case_name;
cfg.n_bodies    = length(body_list);
cfg.water_depth = water_depth;
cfg.fs.r_outer  = r_outer;

domain = struct();
domain.cfg       = cfg;
domain.body_list = body_list;
domain.fs        = mesh_fs;
domain.seabed    = mesh_seabed;
domain.farfield  = mesh_farfield;
domain.waterline = {}; % 无显式 waterline 时留空

%% =========================================================================
%% 4. 直接调用绘图函数
%% =========================================================================
fig = plot_bmf_domain(domain, 'full');