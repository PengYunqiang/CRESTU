function mesh_farfield = generate_farfield_mesh(waterline, cfg, nz_layers)
% GENERATE_FARFIELD_MESH Execute the documented generate_farfield_mesh operation.
%
% Syntax:
%   mesh_farfield = generate_farfield_mesh(waterline, cfg, nz_layers)
%
% Inputs:
%   waterline       : [struct] Ordered waterline nodes in the z = 0 plane, in m.
%   cfg             : [struct] Validated CRESTU configuration with SI-valued physical fields.
%   nz_layers       : [integer scalar or array] Discrete count or index required by the algorithm.
%
% Outputs:
%   mesh_farfield   : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================

    if nargin < 3 || isempty(nz_layers), nz_layers = 8; end

    h = cfg.water_depth;
    if h <= 0, h = cfg.fs.r_outer * 0.8; end

    r_outer   = cfg.fs.r_outer;
    wl_nodes  = waterline.nodes;
    n_pts     = size(wl_nodes, 1);
    is_closed = waterline.is_closed;

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
    if ~is_closed, seg_lens(end) = 0; end

    total_perimeter = sum(seg_lens);
    cum_s = [0; cumsum(seg_lens(1:end-1))];

    if is_closed
        th_start = atan2(wl_nodes(1, 2), wl_nodes(1, 1));
        th_outer = th_start + 2 * pi * (cum_s / total_perimeter);
    else
        th_start = atan2(wl_nodes(1, 2), wl_nodes(1, 1));
        th_end   = atan2(wl_nodes(end, 2), wl_nodes(end, 1));
        th_outer = th_start + (th_end - th_start) * (cum_s / cum_s(end));
    end

    % =====================================================================
    % =====================================================================
    z_lin = linspace(0, -h, nz_layers + 1);

    if is_closed, n_sectors = n_pts; else, n_sectors = n_pts - 1; end
    n_panels = nz_layers * n_sectors;

    vertices = zeros(n_panels, 4, 3);
    centers  = zeros(n_panels, 3);
    normals  = zeros(n_panels, 3);
    areas    = zeros(n_panels, 1);
    e1       = zeros(n_panels, 3);
    e2       = zeros(n_panels, 3);

    p_idx = 0;
    for i = 1:nz_layers
        z_top = z_lin(i);
        z_bot = z_lin(i+1);

        for j = 1:n_sectors
            j_next = j + 1;
            if j_next > n_pts, j_next = 1; end

            th1 = th_outer(j);
            th2 = th_outer(j_next);

            p1 = [r_outer * cos(th1), r_outer * sin(th1), z_top];
            p2 = [r_outer * cos(th2), r_outer * sin(th2), z_top];
            p3 = [r_outer * cos(th2), r_outer * sin(th2), z_bot];
            p4 = [r_outer * cos(th1), r_outer * sin(th1), z_bot];

            p_idx = p_idx + 1;
            vertices(p_idx, 1, :) = p1;
            vertices(p_idx, 2, :) = p2;
            vertices(p_idx, 3, :) = p3;
            vertices(p_idx, 4, :) = p4;

            p_c = (p1 + p2 + p3 + p4) / 4.0;
            centers(p_idx, :) = p_c;

            th_c = (th1 + th2) / 2.0;
            if j_next == 1 && th2 < th1, th_c = (th1 + th2 + 2*pi) / 2.0; end
            normals(p_idx, :) = [-cos(th_c), -sin(th_c), 0.0];

            d13 = p3 - p1; d24 = p4 - p2;
            areas(p_idx) = 0.5 * norm(cross(d13, d24));

            e1(p_idx, :) = [-sin(th_c), cos(th_c), 0.0];
            e2(p_idx, :) = [0.0, 0.0, -1.0];
        end
    end

    mesh_farfield.header      = sprintf('Generic Farfield Wall (R=%.1fm, NP=%d)', r_outer, n_panels);
    mesh_farfield.ulen        = 1.0;
    mesh_farfield.panel_type  = repmat(5, [n_panels, 1]); % 5 = FARFIELD
    mesh_farfield.isx         = cfg.isx;
    mesh_farfield.isy         = cfg.isy;
    mesh_farfield.n_panels    = n_panels;
    mesh_farfield.vertices    = vertices;
    mesh_farfield.centers     = centers;
    mesh_farfield.normals     = normals;
    mesh_farfield.areas       = areas;
    mesh_farfield.e1          = e1;
    mesh_farfield.e2          = e2;
    mesh_farfield.hydrostatics= struct('Vx',0,'Vy',0,'Vz',0,'V_mean',0,'center_of_buoyancy',[0,0,0]);
end
