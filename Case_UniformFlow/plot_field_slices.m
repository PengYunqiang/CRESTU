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