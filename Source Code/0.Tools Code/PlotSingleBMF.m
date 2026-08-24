clear; clc; close all;

%% 1. 文件交互选择
[filename, filepath] = uigetfile({'*.bmf;*.BMF', 'BMF Mesh Files (*.bmf)'; ...
                                  '*.*',         'All Files (*.*)'}, ...
                                 'Select Mesh File');
if isequal(filename, 0)
    disp('User canceled file selection.');
    return;
end

filestr = fullfile(filepath, filename);
fid = fopen(filestr, 'r');
if fid == -1
    error('Cannot open file: %s', filestr);
end

%% 2. 读取网格文件 (支持 BMF / GDF 格式)
[~, ~, ext] = fileparts(filename);

% 读取前 4 行元数据
header_line = fgetl(fid);
line2 = sscanf(fgetl(fid), '%f');
line3 = sscanf(fgetl(fid), '%f');
line4 = sscanf(fgetl(fid), '%f');

if strcmpi(ext, '.bmf')
    ulen = line2(1);
    main_type = 0;
    if numel(line2) >= 2, main_type = line2(2); end
    isx = line3(1); 
    isy = line3(2);
    npan = round(line4(1));
else
    % GDF 格式标准解析 (前4行为报头/重力/对称性等)
    ulen = line2(1);
    npan = round(line4(1));
end

% 读取所有节点坐标数据 (每个 Panel 占 4 行 x 3 列)
data = textscan(fid, '%f');
fclose(fid);

coords = reshape(data{1}, 3, [])';
actual_panels = size(coords, 1) / 4;

if actual_panels ~= npan
    warning('Declared panels (%d) mismatch actual panel count (%d). Using actual count.', npan, actual_panels);
    npan = actual_panels;
end

% 重构为 patch 所需的 (4 x npan) 矩阵
X = reshape(coords(:, 1), 4, npan);
Y = reshape(coords(:, 2), 4, npan);
Z = reshape(coords(:, 3), 4, npan);

%% 3. 计算面元几何中心与单位外法向量
% 面元中心 (4 节点平均)
C_x = mean(X, 1);
C_y = mean(Y, 1);
C_z = mean(Z, 1);

% 向量 1->2 与 1->3 叉乘求法向量 (右手法则)
V12 = [X(2,:) - X(1,:); Y(2,:) - Y(1,:); Z(2,:) - Z(1,:)];
V13 = [X(3,:) - X(1,:); Y(3,:) - Y(1,:); Z(3,:) - Z(1,:)];
N = cross(V12, V13, 1);

% 模长归一化
N_mag = sqrt(sum(N.^2, 1));
N_mag(N_mag < 1e-12) = 1;
N_unit = N ./ N_mag;
Nx = N_unit(1, :);
Ny = N_unit(2, :);
Nz = N_unit(3, :);

%% 4. 可视化渲染
fig = figure('Name', ['Mesh Viewer - ', filename], ...
             'Color', 'w', ...
             'Position', [150, 100, 1050, 800]);

% 4.1 绘制半透明网格面
h_patch = patch(X, Y, Z, [0.3 0.75 0.9], ...
                'FaceAlpha', 0.6, ...
                'EdgeColor', [0.2 0.2 0.2], ...
                'LineWidth', 0.5);

hold on; grid on; axis equal;
view(3);
xlabel('X (m)', 'FontSize', 10);
ylabel('Y (m)', 'FontSize', 10);
zlabel('Z (m)', 'FontSize', 10);

% 4.2 自适应抽样绘制法向量
stride = max(1, round(npan / 50)); % 保持大约显示 50 个法向量
idx = 1:stride:npan;

h_quiver = quiver3(C_x(idx), C_y(idx), C_z(idx), ...
                   Nx(idx), Ny(idx), Nz(idx), ...
                   0.8, 'r', ...
                   'LineWidth', 1.2, ...
                   'Marker', '.', ...
                   'MarkerSize', 8);

% 4.3 标题与图例
title({['Mesh Viewer: ', filename, ' (Total Panels: ', num2str(npan), ')']; ...
       ['\color[rgb]{0.2,0.5,0.8}Mesh Alpha=0.6, \color[rgb]{0.8,0,0}Normals Sampled Every ', ...
       num2str(stride), ' Panels (Scale=0.8)']}, ...
      'Interpreter', 'tex', 'FontSize', 11);

% 4.4 坐标轴显示优化
set(gca, 'Box', 'on', 'LineWidth', 0.8);
rotate3d on;
hold off;
drawnow;