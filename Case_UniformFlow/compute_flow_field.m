function [X_grid, Y_grid, Z_grid, V_field, Cp_field] = compute_flow_field(mesh_body, sigma, U_inf, domain_lim, N_pts)
    np      = mesh_body.n_panels;
    centers = mesh_body.centers;
    normals = mesh_body.normals;
    verts   = mesh_body.vertices;
    e1_all  = mesh_body.e1;
    e2_all  = mesh_body.e2;
    
    x_lin = linspace(-domain_lim, domain_lim, N_pts);
    y_lin = linspace(-domain_lim, domain_lim, N_pts);
    z_lin = linspace(-domain_lim, domain_lim, N_pts);
    [X_grid, Y_grid, Z_grid] = ndgrid(x_lin, y_lin, z_lin);
    
    total_pts = numel(X_grid);
    field_pts = [X_grid(:), Y_grid(:), Z_grid(:)];
    
    V_ind = zeros(total_pts, 3);
    U_mag2 = dot(U_inf, U_inf);
    
    tic;
    parfor i = 1:total_pts
        X_eval = field_pts(i, :);
        v_pt = [0, 0, 0];
        
        for j = 1:np
            P_j = squeeze(verts(j, :, :));
            [u_ij, v_ij, w_ij] = hess_smith_panel_velocity(P_j, X_eval, ...
                                    centers(j, :), e1_all(j, :), e2_all(j, :), normals(j, :));
            v_pt = v_pt + sigma(j) * [u_ij, v_ij, w_ij];
        end
        V_ind(i, :) = v_pt;
    end
    fprintf('>>> 空间点 (%d 点) 速度场计算完成, 耗时: %.2f 秒\n', total_pts, toc);
    
    V_total = repmat(U_inf, total_pts, 1) + V_ind;
    
    % 球体内部点置 NaN 遮罩 (避免内部虚假流场影响可视化)
    R_sph = sqrt(sum(field_pts.^2, 2));
    inside_mask = R_sph < (max(sqrt(sum(centers.^2, 2))) * 0.96);
    V_total(inside_mask, :) = NaN;
    
    Cp_pts = 1.0 - sum(V_total.^2, 2) / U_mag2;
    
    % 重构为 3D 阵列维度
    V_field.u = reshape(V_total(:, 1), size(X_grid));
    V_field.v = reshape(V_total(:, 2), size(X_grid));
    V_field.w = reshape(V_total(:, 3), size(X_grid));
    V_field.mag = reshape(sqrt(sum(V_total.^2, 2)), size(X_grid));
    Cp_field = reshape(Cp_pts, size(X_grid));
end

% =========================================================================
% 2. 经典正交切面可视化函数 (含流线与物面融合)
% =========================================================================
function plot_field_slices(X, Y, Z, V, Cp, mesh_body, diameter)
    fig = figure('Color', 'w', 'Position', [150, 100, 1100, 500], 'Name', 'Flow Field Slices');
    
    % (1) Y = 0 纵剖面对称面流场云图
    subplot(1, 2, 1); hold on; grid on; axis equal;
    mid_y = round(size(Y, 2) / 2);
    x_sub = squeeze(X(:, mid_y, :));
    z_sub = squeeze(Z(:, mid_y, :));
    mag_sub = squeeze(V.mag(:, mid_y, :));
    u_sub = squeeze(V.u(:, mid_y, :));
    w_sub = squeeze(V.w(:, mid_y, :));
    
    pcolor(x_sub, z_sub, mag_sub);
    shading interp;
    colormap(gca, jet);
    c = colorbar;
    c.Label.String = '|V| / U_{\infty}';
    
    % 叠加流线
    [sx, sz] = meshgrid(linspace(min(x_sub(:)), min(x_sub(:)), 15), linspace(min(z_sub(:)), max(z_sub(:)), 15));
    streamline(x_sub', z_sub', u_sub', w_sub', sx, sz);
    
    % 绘制球体截面边界
    th = linspace(0, 2*pi, 100);
    fill(0.5 * diameter * cos(th), 0.5 * diameter * sin(th), [0.8 0.8 0.8], 'EdgeColor', 'k', 'LineWidth', 1.2);
    
    title('Y = 0 纵剖面流速模长云图与流线', 'FontSize', 11, 'FontWeight', 'bold');
    xlabel('X (m)'); ylabel('Z (m)');
    xlim([-diameter, diameter]); ylim([-diameter, diameter]);
    
    % (2) Z = 0 水平赤道剖面流场云图
    subplot(1, 2, 2); hold on; grid on; axis equal;
    mid_z = round(size(Z, 3) / 2);
    x_sub2 = squeeze(X(:, :, mid_z));
    y_sub2 = squeeze(Y(:, :, mid_z));
    mag_sub2 = squeeze(V.mag(:, :, mid_z));
    u_sub2 = squeeze(V.u(:, :, mid_z));
    v_sub2 = squeeze(V.v(:, :, mid_z));
    
    pcolor(x_sub2, y_sub2, mag_sub2);
    shading interp;
    colormap(gca, jet);
    c2 = colorbar;
    c2.Label.String = '|V| / U_{\infty}';
    
    [sx2, sy2] = meshgrid(linspace(min(x_sub2(:)), min(x_sub2(:)), 15), linspace(min(y_sub2(:)), max(y_sub2(:)), 15));
    streamline(x_sub2', y_sub2', u_sub2', v_sub2', sx2, sy2);
    fill(0.5 * diameter * cos(th), 0.5 * diameter * sin(th), [0.8 0.8 0.8], 'EdgeColor', 'k', 'LineWidth', 1.2);
    
    title('Z = 0 水平切面流速模长云图与流线', 'FontSize', 11, 'FontWeight', 'bold');
    xlabel('X (m)'); ylabel('Y (m)');
    xlim([-diameter, diameter]); ylim([-diameter, diameter]);
end