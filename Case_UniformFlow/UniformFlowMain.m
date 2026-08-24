clear; clc; close all;

% =========================================================================
% 1. 生成全封闭球体物面网格 (Full Sphere, D = 10m)
% =========================================================================
mesh_file = 'full_sphere_test.bmf';
diameter = 10.0;
N_grid = 8; % 剖分段数 (6 * 4 * N^2 = 1536 面元)

mesh_body = generate_full_sphere_bmf(mesh_file, diameter, 0, 0, N_grid);

% =========================================================================
% 2. 启动 Hess-Smith 均匀流求解器 (U_inf = 1.0 m/s 沿 +X)
% =========================================================================
U_inf = [1.0, 0.0, 0.0];
flow_res = solve_uniform_flow(mesh_body, U_inf);

% =========================================================================
% 3. 物理守恒与受力对标
% =========================================================================
rho = 1025.0;
p_dyn = 0.5 * rho * dot(U_inf, U_inf) .* flow_res.Cp;
F_hydro = zeros(1, 3);
for i = 1:mesh_body.n_panels
    F_hydro = F_hydro - p_dyn(i) * mesh_body.normals(i, :) * mesh_body.areas(i);
end

fprintf('\n---------------- 均匀流受力与物理守恒对标 ----------------\n');
fprintf(' 设定来流速度 U_inf : [%.2f, %.2f, %.2f] m/s\n', U_inf(1), U_inf(2), U_inf(3));
fprintf(' 表面最大压力系数 Cp_max (驻点理论值 = +1.000) : %8.4f\n', max(flow_res.Cp));
fprintf(' 表面最小压力系数 Cp_min (赤道理论值 = -1.250) : %8.4f\n', min(flow_res.Cp));
fprintf(' 闭合体总阻力 Fx (理论值 = 0 N)                  : %10.4e N\n', F_hydro(1));
fprintf(' 闭合体侧向力 Fy (理论值 = 0 N)                  : %10.4e N\n', F_hydro(2));
fprintf(' 闭合体升力   Fz (理论值 = 0 N)                  : %10.4e N\n', F_hydro(3));
fprintf('========================================================\n');

% =========================================================================
% 4. 提取 Y \approx 0 纵剖线并绘制对标曲线
% =========================================================================
centers = mesh_body.centers;
tol_y = diameter / (2 * N_grid);
idx_slice = find(abs(centers(:, 2)) < tol_y);

x_s = centers(idx_slice, 1);
z_s = centers(idx_slice, 3);
cp_s = flow_res.Cp(idx_slice);

th_deg = atan2d(-z_s, -x_s);
[th_sorted, sort_idx] = sort(th_deg);
cp_sorted = cp_s(sort_idx);

figure('Color', 'w', 'Position', [100, 150, 780, 480], 'Name', 'Cp Validation');
plot(th_sorted, cp_sorted, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 5.5, 'DisplayName', 'BEM 数值解 (Hess-Smith)');
hold on; grid on;

th_theory = linspace(-180, 180, 360);
cp_theory = 1.0 - 2.25 * sind(th_theory).^2;
plot(th_theory, cp_theory, 'b--', 'LineWidth', 2.0, 'DisplayName', '解析解: C_p = 1 - 2.25 sin^2\theta');

xlabel('经向极角 \theta (deg) [0°: 迎流前缘, \pm90°: 赤道, \pm180°: 尾部]', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('表面压力系数 C_p', 'FontSize', 11, 'FontWeight', 'bold');
title('全封闭球体均匀流表面 C_p 对标', 'FontSize', 12);
legend('Location', 'south');
ylim([-1.4, 1.2]);

%% =========================================================================
% 5. 三维流场可视化
% =========================================================================
fig = figure('Color', 'w', 'Position', [900, 150, 850, 650], 'Name', 'Flow 3D');
ax = gca; hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
view(ax, 135, 25);
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('物面压力系数 C_p 云图与流速矢量');

X = squeeze(mesh_body.vertices(:, :, 1))';
Y = squeeze(mesh_body.vertices(:, :, 2))';
Z = squeeze(mesh_body.vertices(:, :, 3))';
colormap(ax, jet);
patch(ax, X, Y, Z, flow_res.Cp', 'FaceColor', 'flat', 'EdgeColor', [0.25 0.25 0.25], 'LineWidth', 0.2);
cbar = colorbar(ax);
cbar.Label.String = 'C_p';

step = 1;
quiver3(ax, mesh_body.centers(1:step:end, 1), ...
            mesh_body.centers(1:step:end, 2), ...
            mesh_body.centers(1:step:end, 3), ...
            flow_res.V_total(1:step:end, 1), ...
            flow_res.V_total(1:step:end, 2), ...
            flow_res.V_total(1:step:end, 3), 1.2, 'k', 'LineWidth', 1.0);

% =========================================================================
% 流场外点诱导速度计算与切片可视化 (调用示例)
% =========================================================================
% 设定计算域为直径的 2 倍范围 ([-D, D])
domain_lim = diameter; 
grid_res   = 41; % 网格分辨率

fprintf('>>> 正在生成流场网格并计算诱导速度场...\n');
[X_grid, Y_grid, Z_grid, V_field, Cp_field] = compute_flow_field(...
    mesh_body, flow_res.sigma, U_inf, domain_lim, grid_res);

% 绘制 Y=0 纵剖面与 Z=0 水平切面的速度模长云图与流线
plot_field_slices(X_grid, Y_grid, Z_grid, V_field, Cp_field, mesh_body, diameter);