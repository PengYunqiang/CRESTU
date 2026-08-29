function mesh = generate_seabed_disk_mesh(cfg, mesh_fs)
% GENERATE_SEABED_DISK_MESH Generate seabed disk mesh for the CRESTU hydrodynamic workflow.
%
% Syntax:
%   mesh = generate_seabed_disk_mesh(cfg, mesh_fs)
%
% Description:
%   Prepares boundary-panel geometry for the Rankine solver.
%   Global coordinates and panel-normal signs are preserved.
%
% Inputs:
%   cfg                - [struct] Validated CRESTU configuration containing SI-valued physical and numerical parameters.
%   mesh_fs            - [struct] Free-surface mesh with coordinates in [m] and areas in [m^2].
%
% Outputs:
%   mesh               - [struct] Generated boundary-panel mesh with coordinates in [m] and areas in [m^2].
%
% Governing Equations / Theory:
%   Planar panel geometry, polygon moments, reflection transformations, and mesh-topology relations as applicable.
%
% References:
%   - Hess and Smith (1964); CRESTU BMF mesh-format and geometry conventions.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    freeSurfaceVertices = reshape(mesh_fs.vertices, mesh_fs.n_panels, 4, 3);
    radialDistance = sqrt(freeSurfaceVertices(:, :, 1) .^ 2 + freeSurfaceVertices(:, :, 2) .^ 2);
    mask = abs(radialDistance - cfg.fs.r_outer) < 1e-6 * max(1, cfg.fs.r_outer);
    xy = [reshape(freeSurfaceVertices(:, :, 1), [], 1), reshape(freeSurfaceVertices(:, :, 2), [], 1)];
    boundary = xy(mask(:), :);
    if size(boundary, 1) < 8
        theta = linspace(0, 2 * pi, 65).';
        theta(end) = [];
        boundary = cfg.fs.r_outer * [cos(theta), sin(theta)];
    else
        angles = atan2(boundary(:, 2), boundary(:, 1));
        [~, order] = sort(angles);
        boundary = unique(boundary(order, :),'rows','stable');
    end
    boundaryPointCount = size(boundary, 1);
    constraints = [(1:boundaryPointCount).', [2:boundaryPointCount, 1].'];
    spacing = max(cfg.fs.r_outer / 14, 0.75);
    a = -cfg.fs.r_outer:spacing:cfg.fs.r_outer;
    [gx, gy] = meshgrid(a, a);
    pts = [gx(:), gy(:)];
    pts = pts(sum(pts .^ 2, 2) < (0.985 * cfg.fs.r_outer)^2, :);
    points = [boundary;pts];
    [points, ~, map] = uniquetol(points, 1e-10,'ByRows', true);
    constraints = map(constraints);
    constraints = constraints(constraints(:, 1) ~= constraints(:, 2), :);
    dt = delaunayTriangulation(points, constraints);
    tri = dt.ConnectivityList;
    points = dt.Points;
    triangleCenters = (points(tri(:, 1), :) + points(tri(:, 2), :) + points(tri(:, 3), :)) / 3;
    tri = tri(sum(triangleCenters .^ 2, 2) <= cfg.fs.r_outer^2, :);
    np = size(tri, 1);
    vertices = zeros(np, 4, 3);
    for vertexIndex = 1:3
        vertices(:, vertexIndex, 1:2) = points(tri(:, vertexIndex), :);
    end
    vertices(:, 4, :) = vertices(:, 3, :);
    vertices(:, :, 3) = -cfg.water_depth;
    [centers, areas, e1, e2] = triangle_geometry(vertices);
    mesh = struct('header', sprintf('Multi-body disk seabed h=%.3g', cfg.water_depth),'ulen', 1, ...
'panel_type', repmat(4, np, 1),'isx', 0,'isy', 0,'n_panels', np,'vertices', vertices, ...
'centers', centers,'normals', repmat([0, 0, 1], np, 1),'areas', areas,'e1', e1,'e2', e2, ...
'hydrostatics', struct('Vx', 0,'Vy', 0,'Vz', 0,'V_mean', 0,'center_of_buoyancy', [0, 0, 0]));
end

function [centers, areas, e1, e2] = triangle_geometry(vertices)
% TRIANGLE_GEOMETRY Evaluate triangular panel centers, areas, and local tangent bases.
%
% Syntax:
%   [centers, areas, e1, e2] = triangle_geometry(vertices)
%
% Description:
%   Prepares boundary-panel geometry for the Rankine solver.
%   Global coordinates and panel-normal signs are preserved.
%
% Inputs:
%   vertices           - [N x 4 x 3] Ordered quadrilateral panel vertices, [m].
%
% Outputs:
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

    np = size(vertices, 1);
    centers = zeros(np, 3);
    areas = zeros(np, 1);
    e1 = zeros(np, 3);
    e2 = zeros(np, 3);
    for panelIndex = 1:np
        p1 = reshape(vertices(panelIndex, 1, :), 1, 3);
        p2 = reshape(vertices(panelIndex, 2, :), 1, 3);
        p3 = reshape(vertices(panelIndex, 3, :), 1, 3);
        centers(panelIndex, :) = (p1 + p2 + p3) / 3;
        areas(panelIndex) = 0.5 * norm(cross(p2 - p1, p3 - p1));
        e1(panelIndex, :) = (p2 - p1) / norm(p2 - p1);
        e2(panelIndex, :) = cross([0, 0, 1], e1(panelIndex, :));
    end
end
