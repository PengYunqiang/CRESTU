clear; clc; close all;

%% --- 1. Select the Input File ---
[filename, filepath] = uigetfile({'*.bmf;*.BMF', 'BMF Mesh Files (*.bmf)'; ...
                                  '*.*', 'All Files (*.*)'}, ...
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

%% --- 2. Read the BMF or GDF Mesh ---
[~, ~, ext] = fileparts(filename);

% Read the first four metadata records.
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
    % Parse the standard GDF header, gravity, and symmetry records.
    ulen = line2(1);
    npan = round(line4(1));
end

% Read all panel vertices as four rows of three coordinates per panel.
data = textscan(fid, '%f');
fclose(fid);

coords = reshape(data{1}, 3, [])';
actual_panels = size(coords, 1) / 4;

if actual_panels ~= npan
    warning('Declared panels (%d) mismatch actual panel count (%d). Using actual count.', npan, actual_panels);
    npan = actual_panels;
end

% Reshape the vertices into 4-by-npan arrays required by patch.
X = reshape(coords(:, 1), 4, npan);
Y = reshape(coords(:, 2), 4, npan);
Z = reshape(coords(:, 3), 4, npan);

%% --- 3. Evaluate Panel Centers and Unit Normals ---
% Evaluate each panel center as the mean of its four vertices.
C_x = mean(X, 1);
C_y = mean(Y, 1);
C_z = mean(Z, 1);

% Evaluate the oriented normal from the 1-to-2 and 1-to-3 edge cross product.
V12 = [X(2, :) - X(1, :); Y(2, :) - Y(1, :); Z(2, :) - Z(1, :)];
V13 = [X(3, :) - X(1, :); Y(3, :) - Y(1, :); Z(3, :) - Z(1, :)];
N = cross(V12, V13, 1);

% Normalize each vector to unit magnitude.
N_mag = sqrt(sum(N .^ 2, 1));
N_mag(N_mag < 1e-12) = 1;
N_unit = N ./ N_mag;
Nx = N_unit(1, :);
Ny = N_unit(2, :);
Nz = N_unit(3, :);

%% --- 4. Render the Mesh ---
fig = figure('Name', ['Mesh Viewer - ', filename], ...
             'Color', 'w', ...
             'Position', [150, 100, 1050, 800]);

% Render translucent panel faces.
h_patch = patch(X, Y, Z, [0.3 0.75 0.9], ...
                'FaceAlpha', 0.6, ...
                'EdgeColor', [0.2 0.2 0.2], ...
                'LineWidth', 0.5);

hold on; grid on; axis equal;
view(3);
xlabel('X (m)', 'FontSize', 10);
ylabel('Y (m)', 'FontSize', 10);
zlabel('Z (m)', 'FontSize', 10);

% Render adaptively sampled panel normals.
stride = max(1, round(npan / 50)); % Retain approximately 50 displayed normal vectors.
idx = 1:stride:npan;

h_quiver = quiver3(C_x(idx), C_y(idx), C_z(idx), ...
                   Nx(idx), Ny(idx), Nz(idx), ...
                   0.8, 'r', ...
                   'LineWidth', 1.2, ...
                   'Marker', '.', ...
                   'MarkerSize', 8);

% Configure the title and legend.
title({['Mesh Viewer: ', filename, ' (Total Panels: ', num2str(npan), ')']; ...
       ['\color[rgb]{0.2,0.5,0.8}Mesh Alpha=0.6, \color[rgb]{0.8,0,0}Normals Sampled Every ', ...
       num2str(stride), ' Panels (Scale=0.8)']}, ...
      'Interpreter', 'tex', 'FontSize', 11);

% Configure equal-axis Cartesian display.
set(gca, 'Box', 'on', 'LineWidth', 0.8);
rotate3d on;
hold off;
drawnow;