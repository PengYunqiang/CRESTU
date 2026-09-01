function mesh_fs = generate_free_surface_bmf(filename, waterline, cfg, writeOutput)
% GENERATE_FREE_SURFACE_BMF Generate free surface bmf for the CRESTU hydrodynamic workflow.
%
% Syntax:
%   mesh_fs = generate_free_surface_bmf(filename, waterline, cfg)
%
% Description:
%   Prepares boundary-panel geometry for the Rankine solver.
%   Global coordinates and panel-normal signs are preserved.
%
% Inputs:
%   filename           - [character vector or string scalar] Input or output file path.
%   waterline          - [struct] Ordered waterline nodes and segment metadata, with coordinates in [m].
%   cfg                - [struct] Validated CRESTU configuration containing SI-valued physical and numerical parameters.
%
% Outputs:
%   mesh_fs            - [struct] Generated free-surface panel mesh in SI units.
%
% Governing Equations / Theory:
%   Planar panel geometry, polygon moments, reflection transformations, and mesh-topology relations as applicable.
%
% References:
%   - Hess and Smith (1964); CRESTU BMF mesh-format and geometry conventions.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    if nargin < 4 || isempty(writeOutput)
        writeOutput = true;
    end
    if nargin < 1 || isempty(filename)
        filename = cfg.files.fs;
    end

    nr_near = cfg.fs.nr_near;
    nr_sponge = cfg.fs.nr_sponge;
    r_inner_base = cfg.fs.r_inner;
    if isfield(cfg, 'phase2_2_controls') && ...
            cfg.phase2_2_controls.sectionPresent
        r_inner_base = cfg.phase2_2_controls.meshTransitionRadiusM;
    end
    r_outer = cfg.fs.r_outer;
    sponge_ratio = cfg.fs.sponge_ratio;
    mu0 = cfg.fs.mu0;

    wl_nodes = waterline.nodes; % [N_pts x 2]
    n_pts = size(wl_nodes, 1);
    is_closed = waterline.is_closed;
    nr_total = nr_near + nr_sponge;

    if get_polar_topology_mode(cfg) == "QUALITY_CONTROLLED_V2"
        mesh_fs = generate_quality_controlled_polar_surface( ...
            waterline, cfg, 0.0, 2);
        if writeOutput
            write_bmf(filename, mesh_fs);
        end
        return
    end

% =====================================================================
% =====================================================================
    seg_lens = zeros(n_pts, 1);
    for j = 1:n_pts
        if is_closed
            j_next = mod(j, n_pts) + 1;
        else
            j_next = min(j + 1, n_pts);
        end
        seg_lens(j) = norm(wl_nodes(j_next, :) - wl_nodes(j, :));
    end
    if ~is_closed
        seg_lens(end) = 0;
    end

    total_perimeter = sum(seg_lens);
    cum_s = [0; cumsum(seg_lens(1:end - 1))];

    if is_closed
        th_start = atan2(wl_nodes(1, 2), wl_nodes(1, 1));
        th_outer = th_start + 2 * pi * (cum_s / total_perimeter);
    else
        th_start = atan2(wl_nodes(1, 2), wl_nodes(1, 1));
        th_end = atan2(wl_nodes(end, 2), wl_nodes(end, 1));
        th_outer = th_start + (th_end - th_start) * (cum_s / cum_s(end));
    end

    outer_pts = [r_outer * cos(th_outer), r_outer * sin(th_outer)];

% =====================================================================
% =====================================================================
% The old implementation hard-coded the near region to 55% of the
% entire radius.  Once Rout was enlarged, damping began inside an
% under-resolved near field.  Build each ray so row NR_near+1 lies
% exactly on R_inner, followed by a geometrically stretched sponge.
    near_coordinate = linspace(0, 1, nr_near + 1);
    sp_ratios = sponge_ratio .^ (0:nr_sponge - 1);
    sponge_coordinate = cumsum(sp_ratios / sum(sp_ratios));

    grid_x = zeros(nr_total + 1, n_pts);
    grid_y = zeros(nr_total + 1, n_pts);

    for j = 1:n_pts
        local_radius = norm(wl_nodes(j, :));
        interface_fraction = (r_inner_base - local_radius) / (r_outer - local_radius);
        interface_fraction = max(0, min(1, interface_fraction));
        weights = [interface_fraction * near_coordinate, ...
            interface_fraction + (1 - interface_fraction) * sponge_coordinate];
        grid_x(:, j) = (1 - weights(:)) * wl_nodes(j, 1) + weights(:) * outer_pts(j, 1);
        grid_y(:, j) = (1 - weights(:)) * wl_nodes(j, 2) + weights(:) * outer_pts(j, 2);
    end

% =====================================================================
% =====================================================================
    max_iter = 250;
    omega = 1.15;

    for iter = 1:max_iter
        for i = 2:nr_total
            if i == nr_near + 1
                continue;
            end
            for j = 1:n_pts
                if is_closed
                    j_m = mod(j - 2 + n_pts, n_pts) + 1;
                    j_p = mod(j, n_pts) + 1;
                else
                    if j == 1
                        j_m = 2;
                        j_p = 2;
                    elseif j == n_pts
                        j_m = n_pts - 1;
                        j_p = n_pts - 1;
                    else
                        j_m = j - 1;
                        j_p = j + 1;
                    end
                end

                x_xi = (grid_x(i + 1, j) - grid_x(i - 1, j)) * 0.5;
                y_xi = (grid_y(i + 1, j) - grid_y(i - 1, j)) * 0.5;
                x_eta = (grid_x(i, j_p) - grid_x(i, j_m)) * 0.5;
                y_eta = (grid_y(i, j_p) - grid_y(i, j_m)) * 0.5;

                alpha = x_eta^2 + y_eta^2;
                gamma = x_xi^2 + y_xi^2;
                beta = x_xi * x_eta + y_xi * y_eta;

                if is_closed
                    j_pp = mod(j, n_pts) + 1;
                    j_mm = mod(j - 2 + n_pts, n_pts) + 1;
                else
                    j_pp = min(j + 1, n_pts);
                    j_mm = max(j - 1, 1);
                end
                x_xi_eta = (grid_x(i + 1, j_pp) - grid_x(i - 1, j_pp) - grid_x(i + 1, j_mm) + grid_x(i - 1, j_mm)) * 0.25;
                y_xi_eta = (grid_y(i + 1, j_pp) - grid_y(i - 1, j_pp) - grid_y(i + 1, j_mm) + grid_y(i - 1, j_mm)) * 0.25;

                denom = 2 * (alpha + gamma + 1e-10);

                x_target = (alpha * (grid_x(i + 1, j) + grid_x(i - 1, j)) + ...
                            gamma * (grid_x(i, j_p) + grid_x(i, j_m)) - ...
                            2 * beta * x_xi_eta) / denom;

                y_target = (alpha * (grid_y(i + 1, j) + grid_y(i - 1, j)) + ...
                            gamma * (grid_y(i, j_p) + grid_y(i, j_m)) - ...
                            2 * beta * y_xi_eta) / denom;

                grid_x(i, j) = (1 - omega) * grid_x(i, j) + omega * x_target;
                grid_y(i, j) = (1 - omega) * grid_y(i, j) + omega * y_target;
            end
        end
    end

% =====================================================================
% =====================================================================
    if is_closed
        n_sectors = n_pts;
    else
        n_sectors = n_pts - 1;
    end
    n_panels = nr_total * n_sectors;

    vertices = zeros(n_panels, 4, 3);
    centers = zeros(n_panels, 3);
    normals = zeros(n_panels, 3);
    areas = zeros(n_panels, 1);
    mu_layer = zeros(n_panels, 1);
    e1 = zeros(n_panels, 3);
    e2 = zeros(n_panels, 3);

    p_idx = 0;
    for i = 1:nr_total
        for j = 1:n_sectors
            j_next = j + 1;
            if j_next > n_pts
                j_next = 1;
            end

            p1 = [grid_x(i, j), grid_y(i, j), 0.0];
            p2 = [grid_x(i + 1, j), grid_y(i + 1, j), 0.0];
            p3 = [grid_x(i + 1, j_next), grid_y(i + 1, j_next), 0.0];
            p4 = [grid_x(i, j_next), grid_y(i, j_next), 0.0];

            p_idx = p_idx + 1;
            vertices(p_idx, 1, :) = p1;
            vertices(p_idx, 2, :) = p2;
            vertices(p_idx, 3, :) = p3;
            vertices(p_idx, 4, :) = p4;

            p_c = (p1 + p2 + p3 + p4) / 4.0;
            centers(p_idx, :) = p_c;
            normals(p_idx, :) = [0.0, 0.0, 1.0];

            d13 = p3 - p1;
            d24 = p4 - p2;
            areas(p_idx) = 0.5 * norm(cross(d13, d24));

            v_rad = (p2 + p3) / 2 - (p1 + p4) / 2;
            e1_vec = v_rad / norm(v_rad);
            e1(p_idx, :) = e1_vec;
            e2(p_idx, :) = cross([0, 0, 1], e1_vec);

            r_c = norm(p_c(1:2));
            if r_c > r_inner_base
                mu_layer(p_idx) = mu0 * ((r_c - r_inner_base) / (r_outer - r_inner_base))^2;
            else
                mu_layer(p_idx) = 0.0;
            end
        end
    end

    mesh_fs.header = sprintf('Elliptic Smoothed FS (NP=%d, Rout=%.1fm)', n_panels, r_outer);
    mesh_fs.ulen = 1.0;
    mesh_fs.panel_type = repmat(2, [n_panels, 1]);
    mesh_fs.isx = cfg.isx;
    mesh_fs.isy = cfg.isy;
    mesh_fs.n_panels = n_panels;
    mesh_fs.vertices = vertices;
    mesh_fs.centers = centers;
    mesh_fs.normals = normals;
    mesh_fs.areas = areas;
    mesh_fs.e1 = e1;
    mesh_fs.e2 = e2;
    mesh_fs.mu_damping = mu_layer;
    mesh_fs.hydrostatics = struct('Vx', 0,'Vy', 0,'Vz', 0,'V_mean', 0,'center_of_buoyancy', [0, 0, 0]);

    if writeOutput
        write_bmf(filename, mesh_fs);
    end
end
