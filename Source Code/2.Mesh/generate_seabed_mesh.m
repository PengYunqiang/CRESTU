function mesh_seabed = generate_seabed_mesh(waterline, cfg, mesh_fs)
% GENERATE_SEABED_MESH Generate seabed mesh for the CRESTU hydrodynamic workflow.
%
% Syntax:
%   mesh_seabed = generate_seabed_mesh(waterline, cfg, mesh_fs)
%
% Description:
%   Prepares boundary-panel geometry for the Rankine solver.
%   Global coordinates and panel-normal signs are preserved.
%
% Inputs:
%   waterline          - [struct] Ordered waterline nodes and segment metadata, with coordinates in [m].
%   cfg                - [struct] Validated CRESTU configuration containing SI-valued physical and numerical parameters.
%   mesh_fs            - [struct] Free-surface mesh with coordinates in [m] and areas in [m^2].
%
% Outputs:
%   mesh_seabed        - [struct] Generated seabed panel mesh in SI units.
%
% Governing Equations / Theory:
%   Planar panel geometry, polygon moments, reflection transformations, and mesh-topology relations as applicable.
%
% References:
%   - Hess and Smith (1964); CRESTU BMF mesh-format and geometry conventions.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    waterDepth = cfg.water_depth;
    if waterDepth <= 0
        error('A seabed mesh is not required for an infinite-depth configuration');
    end

    if nargin < 3 || isempty(mesh_fs)
        mesh_fs = generate_free_surface_bmf('', waterline, cfg);
    end

% =====================================================================
% =====================================================================
    n_ring_panels = mesh_fs.n_panels;
    v_ring = mesh_fs.vertices;
    v_ring(:, :, 3) = -waterDepth;
    c_ring = mesh_fs.centers;
    c_ring(:, 3) = -waterDepth;
    a_ring = mesh_fs.areas;
    e1_ring = mesh_fs.e1;
    e2_ring = mesh_fs.e2;

% =====================================================================
% =====================================================================
    wl_nodes = waterline.nodes;
    n_pts = size(wl_nodes, 1);
    wl_closed = [wl_nodes; wl_nodes(1, :)];

    xc = mean(wl_nodes(:, 1));
    yc = mean(wl_nodes(:, 2));

    d_vec = diff(wl_closed, 1, 1);
    d_len = sqrt(sum(d_vec .^ 2, 2)) + 1e-12;
    d_norm = d_vec ./ repmat(d_len, 1, 2);

    turn_angles = zeros(n_pts, 1);
    for i = 1:n_pts
        i_prev = mod(i - 2 + n_pts, n_pts) + 1;
        cos_th = dot(d_norm(i_prev, :), d_norm(i, :));
        cos_th = max(-1.0, min(1.0, cos_th));
        turn_angles(i) = acos(cos_th);
    end

    [sorted_angles, sort_idx] = sort(turn_angles,'descend');

    if sorted_angles(4) > deg2rad(30)
        c_idx = sort(sort_idx(1:4));
    else
        c_idx = round(linspace(1, n_pts + 1, 5));
        c_idx = c_idx(1:4)';
    end

    idx_next = [c_idx(2:4); c_idx(1)];
    seg_counts = zeros(4, 1);
    for k = 1:4
        if idx_next(k) > c_idx(k)
            seg_counts(k) = idx_next(k) - c_idx(k);
        else
            seg_counts(k) = idx_next(k) + n_pts - c_idx(k);
        end
    end

    Nx = max(2, round((seg_counts(1) + seg_counts(3)) / 2));
    Ny = max(2, round((seg_counts(2) + seg_counts(4)) / 2));

    seg1_pts = extract_resample_segment(wl_nodes, c_idx(1), c_idx(2), Nx + 1);
    seg2_pts = extract_resample_segment(wl_nodes, c_idx(2), c_idx(3), Ny + 1);
    seg3_pts = extract_resample_segment(wl_nodes, c_idx(3), c_idx(4), Nx + 1);
    seg4_pts = extract_resample_segment(wl_nodes, c_idx(4), c_idx(1), Ny + 1);

% =====================================================================
% =====================================================================
    C1 = seg1_pts(1, :);
    C2 = seg2_pts(1, :);
    C3 = seg3_pts(1, :);
    C4 = seg4_pts(1, :);

    diskCenter = [xc, yc];
    alpha = 0.5;

    K1 = diskCenter + alpha * (C1 - diskCenter);
    K2 = diskCenter + alpha * (C2 - diskCenter);
    K3 = diskCenter + alpha * (C3 - diskCenter);
    K4 = diskCenter + alpha * (C4 - diskCenter);

    edge1 = [linspace(K1(1), K2(1), Nx + 1)', linspace(K1(2), K2(2), Nx + 1)'];
    edge2 = [linspace(K2(1), K3(1), Ny + 1)', linspace(K2(2), K3(2), Ny + 1)'];
    edge3 = [linspace(K3(1), K4(1), Nx + 1)', linspace(K3(2), K4(2), Nx + 1)'];
    edge4 = [linspace(K4(1), K1(1), Ny + 1)', linspace(K4(2), K1(2), Ny + 1)'];

% =====================================================================
% =====================================================================
    Nr = max(3, round((Nx + Ny) * 0.25));
    requestedQualityNr = 0;
    if get_polar_topology_mode(cfg) == "QUALITY_CONTROLLED_V2"
        assert(isfield(cfg.phase2_2_controls, ...
            'bottomCoreRadialLayerCount'), ...
            'CRESTU:QualityBottomCoreControlMissing', ...
            'QUALITY_CONTROLLED_V2 requires BOTTOM_CORE_RADIAL_LAYERS.');
        requestedQualityNr = ...
            cfg.phase2_2_controls.bottomCoreRadialLayerCount;
        assert(requestedQualityNr >= 1 && fix(requestedQualityNr) == ...
            requestedQualityNr, 'CRESTU:QualityBottomCoreControl', ...
            'BOTTOM_CORE_RADIAL_LAYERS must be a positive integer.');
        Nr = requestedQualityNr;
    end
    radialBlend = linspace(0.0, 1.0, Nr + 1);

% Block 1
    G1_x = zeros(Nx + 1, Nr + 1);
    G1_y = zeros(Nx + 1, Nr + 1);
    for k = 1:Nr + 1
        G1_x(:, k) = (1 - radialBlend(k)) * seg1_pts(:, 1) + radialBlend(k) * edge1(:, 1);
        G1_y(:, k) = (1 - radialBlend(k)) * seg1_pts(:, 2) + radialBlend(k) * edge1(:, 2);
    end

% Block 2
    G2_x = zeros(Ny + 1, Nr + 1);
    G2_y = zeros(Ny + 1, Nr + 1);
    for k = 1:Nr + 1
        G2_x(:, k) = (1 - radialBlend(k)) * seg2_pts(:, 1) + radialBlend(k) * edge2(:, 1);
        G2_y(:, k) = (1 - radialBlend(k)) * seg2_pts(:, 2) + radialBlend(k) * edge2(:, 2);
    end

% Block 3
    G3_x = zeros(Nx + 1, Nr + 1);
    G3_y = zeros(Nx + 1, Nr + 1);
    for k = 1:Nr + 1
        G3_x(:, k) = (1 - radialBlend(k)) * seg3_pts(:, 1) + radialBlend(k) * edge3(:, 1);
        G3_y(:, k) = (1 - radialBlend(k)) * seg3_pts(:, 2) + radialBlend(k) * edge3(:, 2);
    end

% Block 4
    G4_x = zeros(Ny + 1, Nr + 1);
    G4_y = zeros(Ny + 1, Nr + 1);
    for k = 1:Nr + 1
        G4_x(:, k) = (1 - radialBlend(k)) * seg4_pts(:, 1) + radialBlend(k) * edge4(:, 1);
        G4_y(:, k) = (1 - radialBlend(k)) * seg4_pts(:, 2) + radialBlend(k) * edge4(:, 2);
    end

% =====================================================================
% =====================================================================
    u_lin = linspace(0.0, 1.0, Nx + 1);
    v_lin = linspace(0.0, 1.0, Ny + 1);
    [U, V] = ndgrid(u_lin, v_lin);

    G0_x = (1 - U) .* (1 - V) * K1(1) + U .* (1 - V) * K2(1) + ...
            U .* V * K3(1) + (1 - U) .* V * K4(1);
    G0_y = (1 - U) .* (1 - V) * K1(2) + U .* (1 - V) * K2(2) + ...
            U .* V * K3(2) + (1 - U) .* V * K4(2);

% =====================================================================
% =====================================================================
    G0_x = smooth_internal_grid(G0_x);
    G0_y = smooth_internal_grid(G0_y);
    G1_x = smooth_internal_grid(G1_x);
    G1_y = smooth_internal_grid(G1_y);
    G2_x = smooth_internal_grid(G2_x);
    G2_y = smooth_internal_grid(G2_y);
    G3_x = smooth_internal_grid(G3_x);
    G3_y = smooth_internal_grid(G3_y);
    G4_x = smooth_internal_grid(G4_x);
    G4_y = smooth_internal_grid(G4_y);

% =====================================================================
% =====================================================================
    [v1, c1, a1, e1_1, e2_1] = grid_to_panels(G1_x, G1_y, -waterDepth);
    [v2, c2, a2, e1_2, e2_2] = grid_to_panels(G2_x, G2_y, -waterDepth);
    [v3, c3, a3, e1_3, e2_3] = grid_to_panels(G3_x, G3_y, -waterDepth);
    [v4, c4, a4, e1_4, e2_4] = grid_to_panels(G4_x, G4_y, -waterDepth);
    [v0, c0, a0, e1_0, e2_0] = grid_to_panels(G0_x, G0_y, -waterDepth);

    v_infill = cat(1, v1, v2, v3, v4, v0);
    c_infill = cat(1, c1, c2, c3, c4, c0);
    a_infill = cat(1, a1, a2, a3, a4, a0);
    e1_infill = cat(1, e1_1, e1_2, e1_3, e1_4, e1_0);
    e2_infill = cat(1, e2_1, e2_2, e2_3, e2_4, e2_0);

    total_panels = n_ring_panels + size(v_infill, 1);
    vertices = cat(1, v_ring, v_infill);
    centers = cat(1, c_ring, c_infill);
    areas = cat(1, a_ring, a_infill);
    e1 = cat(1, e1_ring, e1_infill);
    e2 = cat(1, e2_ring, e2_infill);
    normals = repmat([0.0, 0.0, 1.0], [total_panels, 1]);

    mesh_seabed.header = sprintf('Adaptive 5-Block Seabed (Depth=%.2fm, NP=%d, Nx=%d, Ny=%d)', ...
                                      waterDepth, total_panels, Nx, Ny);
    mesh_seabed.ulen = 1.0;
    mesh_seabed.panel_type = repmat(4, [total_panels, 1]);
    mesh_seabed.isx = cfg.isx;
    mesh_seabed.isy = cfg.isy;
    mesh_seabed.n_panels = total_panels;
    mesh_seabed.vertices = vertices;
    mesh_seabed.centers = centers;
    mesh_seabed.normals = normals;
    mesh_seabed.areas = areas;
    mesh_seabed.e1 = e1;
    mesh_seabed.e2 = e2;
    mesh_seabed.hydrostatics = struct('Vx', 0,'Vy', 0,'Vz', 0,'V_mean', 0,'center_of_buoyancy', [0, 0, 0]);
    mesh_seabed.quality_controlled_bottom_core = struct( ...
        'topologyMode', char(get_polar_topology_mode(cfg)), ...
        'boundarySegmentCountX', Nx, 'boundarySegmentCountY', Ny, ...
        'legacyRadialLayerCount', max(3, round((Nx + Ny) * 0.25)), ...
        'requestedRadialLayerCount', requestedQualityNr, ...
        'realizedRadialLayerCount', Nr, ...
        'sideBlockPanelCount', 2 * Nx * Nr + 2 * Ny * Nr, ...
        'centralPatchPanelCount', Nx * Ny, ...
        'totalStructuredCorePanelCount', ...
            2 * Nx * Nr + 2 * Ny * Nr + Nx * Ny, ...
        'singleCenterFan', false, 'boundaryNodeContractChanged', false);
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function pts_out = extract_resample_segment(wl_nodes, idx_start, idx_end, n_target)
% EXTRACT_RESAMPLE_SEGMENT Extract and uniformly resample an ordered waterline segment.
%
% Syntax:
%   pts_out = extract_resample_segment(wl_nodes, idx_start, idx_end, n_target)
%
% Description:
%   Prepares boundary-panel geometry for the Rankine solver.
%   Global coordinates and panel-normal signs are preserved.
%
% Inputs:
%   wl_nodes           - [N x 2] Ordered waterline nodes, [m].
%   idx_start          - [scalar] One-based starting node index, dimensionless.
%   idx_end            - [scalar] One-based ending node index, dimensionless.
%   n_target           - [scalar] Requested number of resampled points, dimensionless.
%
% Outputs:
%   pts_out            - [Ntarget x 2] Uniformly resampled planar segment points, [m].
%
% Governing Equations / Theory:
%   Planar panel geometry, polygon moments, reflection transformations, and mesh-topology relations as applicable.
%
% References:
%   - Hess and Smith (1964); CRESTU BMF mesh-format and geometry conventions.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    n_pts = size(wl_nodes, 1);

    if idx_end > idx_start
        indices = idx_start:idx_end;
    else
        indices = [idx_start:n_pts, 1:idx_end];
    end

    raw_pts = wl_nodes(indices, :);
    seg_lens = sqrt(sum(diff(raw_pts, 1, 1) .^ 2, 2));
    cum_s = [0; cumsum(seg_lens)];
    total_len = cum_s(end);

    if total_len < 1e-8
        pts_out = repmat(raw_pts(1, :), n_target, 1);
        return;
    end

    s_target = linspace(0, total_len, n_target);
    pts_out = [interp1(cum_s, raw_pts(:, 1), s_target,'linear')', ...
               interp1(cum_s, raw_pts(:, 2), s_target,'linear')'];
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function G = smooth_internal_grid(G)
% SMOOTH_INTERNAL_GRID Smooth the internal nodes of a structured surface grid.
%
% Syntax:
%   G = smooth_internal_grid(G)
%
% Description:
%   Prepares boundary-panel geometry for the Rankine solver.
%   Global coordinates and panel-normal signs are preserved.
%
% Inputs:
%   G                  - [M x N x 2] Smoothed structured planar grid coordinates, [m].
%
% Outputs:
%   G                  - [M x N x 2] Smoothed structured planar grid coordinates, [m].
%
% Governing Equations / Theory:
%   Planar panel geometry, polygon moments, reflection transformations, and mesh-topology relations as applicable.
%
% References:
%   - Hess and Smith (1964); CRESTU BMF mesh-format and geometry conventions.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    [Ni, Nj] = size(G);
    if Ni <= 2 || Nj <= 2
        return;
    end
    for iter = 1:25
        for i = 2:Ni - 1
            for j = 2:Nj - 1
                G(i, j) = 0.25 * (G(i + 1, j) + G(i - 1, j) + G(i, j + 1) + G(i, j - 1));
            end
        end
    end
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function [vertices, centers, areas, e1, e2] = grid_to_panels(Gx, Gy, z_val)
% GRID_TO_PANELS Convert a structured node grid into quadrilateral panel geometry.
%
% Syntax:
%   [vertices, centers, areas, e1, e2] = grid_to_panels(Gx, Gy, z_val)
%
% Description:
%   Prepares boundary-panel geometry for the Rankine solver.
%   Global coordinates and panel-normal signs are preserved.
%
% Inputs:
%   Gx                 - [M x N] Structured-grid x coordinates, [m].
%   Gy                 - [M x N] Structured-grid y coordinates, [m].
%   z_val              - [scalar] Constant panel-grid elevation, [m].
%
% Outputs:
%   vertices           - [N x 4 x 3] Ordered quadrilateral panel vertices, [m].
%   centers            - [N x 3] Panel collocation points in global coordinates, [m].
%   areas              - [N x 1] Panel areas, [m^2].
%   e1                 - [1 x 3] First unit tangent vector, dimensionless.
%   e2                 - [1 x 3] Second unit tangent vector, dimensionless.
%
% Governing Equations / Theory:
%   Planar panel geometry, polygon moments, reflection transformations, and mesh-topology relations as applicable.
%
% References:
%   - Hess and Smith (1964); CRESTU BMF mesh-format and geometry conventions.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    [Ni, Nj] = size(Gx);
    n_panels = (Ni - 1) * (Nj - 1);

    vertices = zeros(n_panels, 4, 3);
    centers = zeros(n_panels, 3);
    areas = zeros(n_panels, 1);
    e1 = zeros(n_panels, 3);
    e2 = zeros(n_panels, 3);

    p_idx = 0;
    for i = 1:Ni - 1
        for j = 1:Nj - 1
            p_idx = p_idx + 1;

            p1 = [Gx(i, j), Gy(i, j), z_val];
            p2 = [Gx(i + 1, j), Gy(i + 1, j), z_val];
            p3 = [Gx(i + 1, j + 1), Gy(i + 1, j + 1), z_val];
            p4 = [Gx(i, j + 1), Gy(i, j + 1), z_val];

            vertices(p_idx, 1, :) = p1;
            vertices(p_idx, 2, :) = p2;
            vertices(p_idx, 3, :) = p3;
            vertices(p_idx, 4, :) = p4;

            centers(p_idx, :) = (p1 + p2 + p3 + p4) / 4.0;
            d13 = p3 - p1;
            d24 = p4 - p2;
            areas(p_idx) = 0.5 * norm(cross(d13, d24));

            e1(p_idx, :) = [1.0, 0.0, 0.0];
            e2(p_idx, :) = [0.0, 1.0, 0.0];
        end
    end
end
