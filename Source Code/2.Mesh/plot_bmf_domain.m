function fig = plot_bmf_domain(domain, mode)
% PLOT_BMF_DOMAIN Plot bmf domain for the CRESTU hydrodynamic workflow.
%
% Syntax:
%   fig = plot_bmf_domain(domain, mode)
%
% Description:
%   The routine constructs, transforms, validates, or visualizes boundary-panel geometry used by the Rankine solver. Coordinates are expressed in the global Cartesian frame and panel orientation is preserved so that normals remain consistent with boundary-integral signs.
%
% Inputs:
%   domain             - [struct] Assembled body, free-surface, seabed, and far-field boundary domain in SI units.
%   mode               - [character vector or string scalar] Domain-plot rendering mode.
%
% Outputs:
%   fig                - [graphics handle] MATLAB figure containing the requested visualization.
%
% Governing Equations / Theory:
%   Planar panel geometry, polygon moments, reflection transformations, and mesh-topology relations as applicable.
%
% References:
%   - Hess and Smith (1964); CRESTU BMF mesh-format and geometry conventions.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% --- 1. Validate Inputs and Initialize the Algorithm ---

    if nargin < 2 || isempty(mode)
        mode = 'full';
    end

    cfg = domain.cfg;
    body_list = domain.body_list;
    mesh_fs = domain.fs;
    waterline = domain.waterline;
    has_seabed = (cfg.water_depth > 0);

    fig = figure('Color', 'w', 'Position', [100, 100, 1050, 700], ...
                 'Name', sprintf('Rankine Domain [%s] - Mode: %s', cfg.case_name, mode));
    ax = gca;
    hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
    view(ax, 138, 28);
    xlabel(ax, 'X (m)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel(ax, 'Y (m)', 'FontSize', 11, 'FontWeight', 'bold');
    zlabel(ax, 'Z (m)', 'FontSize', 11, 'FontWeight', 'bold');
    title(ax, sprintf('Rankine BEM Fluid Domain: [%s] (Mode: %s)', cfg.case_name, mode), ...
          'FontSize', 13, 'Interpreter', 'none');

    is_wireframe = strcmpi(mode, 'wireframe');
    show_normals = strcmpi(mode, 'full');

    for b = 1:cfg.n_bodies
        mb = body_list{b};
        [Xb, Yb, Zb] = extract_patch_coords(mb);

        if is_wireframe
            patch(ax, Xb, Yb, Zb, 'w', 'FaceColor', 'none', ...
                  'EdgeColor', [0.1 0.25 0.7], 'LineWidth', 0.6);
        else
            patch(ax, Xb, Yb, Zb, 'r', 'FaceColor', [0.25 0.55 0.85], ...
                  'EdgeColor', [0.1 0.25 0.5], 'LineWidth', 0.4, 'FaceAlpha', 0.9);
        end

        if show_normals
            step_b = max(1, round(mb.n_panels / 25));
            quiver3(ax, mb.centers(1:step_b:end, 1), mb.centers(1:step_b:end, 2), mb.centers(1:step_b:end, 3), ...
                        mb.normals(1:step_b:end, 1), mb.normals(1:step_b:end, 2), mb.normals(1:step_b:end, 3), ...
                        1.2, 'r', 'LineWidth', 1.2, 'MaxHeadSize', 0.8);
        end
    end

    [Xf, Yf, Zf] = extract_patch_coords(mesh_fs);
    if is_wireframe
        patch(ax, Xf, Yf, Zf, 'w', 'FaceColor', 'none', ...
              'EdgeColor', [0.3 0.6 0.6], 'LineWidth', 0.3);
    else
        colormap(ax, parula);
        C_fs = reshape(mesh_fs.mu_damping, [1, mesh_fs.n_panels]);
        patch(ax, Xf, Yf, Zf, C_fs, 'CDataMapping', 'scaled', ...
              'EdgeColor', [0.35 0.45 0.45], 'LineWidth', 0.3, 'FaceAlpha', 0.65);

        cbar = colorbar(ax, 'Location', 'eastoutside');
        cbar.Label.String = 'Rayleigh Damping \mu(r) (Sponge Layer)';
        cbar.Label.FontSize = 10;
    end

    if has_seabed
        [Xs, Ys, Zs] = extract_patch_coords(domain.seabed);
        if is_wireframe
            patch(ax, Xs, Ys, Zs, 'w', 'FaceColor', 'none', ...
                  'EdgeColor', [0.6 0.4 0.2], 'LineWidth', 0.3);
        else
            patch(ax, Xs, Ys, Zs, 'r', 'FaceColor', [0.72 0.55 0.35], ...
                  'EdgeColor', [0.35 0.25 0.15], 'LineWidth', 0.3, 'FaceAlpha', 0.5);
        end

        [Xw, Yw, Zw] = extract_patch_coords(domain.farfield);
        if is_wireframe
            patch(ax, Xw, Yw, Zw, 'w', 'FaceColor', 'none', ...
                  'EdgeColor', [0.6 0.3 0.6], 'LineWidth', 0.3);
        else
            patch(ax, Xw, Yw, Zw, 'r', 'FaceColor', [0.75 0.45 0.75], ...
                  'EdgeColor', [0.4 0.2 0.4], 'LineWidth', 0.3, 'FaceAlpha', 0.35);
        end
    end

    if ~iscell(waterline), waterline = {waterline}; end
    for b = 1:numel(waterline)
        wl = waterline{b};
        if isempty(wl) || ~isfield(wl, 'nodes') || isempty(wl.nodes), continue; end
        if wl.is_closed
            wl_pts = zeros(size(wl.nodes, 1) + 1, 2);
            wl_pts(1:end - 1, :) = wl.nodes;
            wl_pts(end, :) = wl.nodes(1, :);
        else
            wl_pts = wl.nodes;
        end
        plot3(ax, wl_pts(:, 1), wl_pts(:, 2), zeros(size(wl_pts, 1), 1), ...
            'r-o', 'LineWidth', 2.2, 'MarkerSize', 3.5, 'MarkerFaceColor', 'r');
    end

    xlim(ax, [-cfg.fs.r_outer * 1.05, cfg.fs.r_outer * 1.05]);
    ylim(ax, [-cfg.fs.r_outer * 1.05, cfg.fs.r_outer * 1.05]);
    if has_seabed
        zlim(ax, [-cfg.water_depth * 1.1, 2.0]);
    else
        zlim(ax, [-10.0, 2.0]);
    end
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function [X, Y, Z] = extract_patch_coords(m)
% EXTRACT_PATCH_COORDS Convert panel vertices to coordinate arrays for patch rendering.
%
% Syntax:
%   [X, Y, Z] = extract_patch_coords(m)
%
% Description:
%   The routine constructs, transforms, validates, or visualizes boundary-panel geometry used by the Rankine solver. Coordinates are expressed in the global Cartesian frame and panel orientation is preserved so that normals remain consistent with boundary-integral signs.
%
% Inputs:
%   m                  - [struct] Boundary mesh containing panel vertices, centers, normals, and areas in SI units.
%
% Outputs:
%   X                  - [4 x N] Panel x coordinates for patch rendering, [m].
%   Y                  - [4 x N] Panel y coordinates for patch rendering, [m].
%   Z                  - [4 x N] Panel z coordinates for patch rendering, [m].
%
% Governing Equations / Theory:
%   Planar panel geometry, polygon moments, reflection transformations, and mesh-topology relations as applicable.
%
% References:
%   - Hess and Smith (1964); CRESTU BMF mesh-format and geometry conventions.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% --- 1. Validate Inputs and Initialize the Algorithm ---

    np = m.n_panels;
    X = squeeze(m.vertices(:, :, 1))';
    Y = squeeze(m.vertices(:, :, 2))';
    Z = squeeze(m.vertices(:, :, 3))';

    if np == 1
        X = X(:); Y = Y(:); Z = Z(:);
    end
end
