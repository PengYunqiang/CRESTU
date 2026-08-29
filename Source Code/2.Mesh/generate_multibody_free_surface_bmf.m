function mesh = generate_multibody_free_surface_bmf(filename, waterlines, cfg)
% GENERATE_MULTIBODY_FREE_SURFACE_BMF Generate multibody free surface bmf for the CRESTU hydrodynamic workflow.
%
% Syntax:
%   mesh = generate_multibody_free_surface_bmf(filename, waterlines, cfg)
%
% Description:
%   Prepares boundary-panel geometry for the Rankine solver.
%   Global coordinates and panel-normal signs are preserved.
%
% Inputs:
%   filename           - [character vector or string scalar] Input or output file path.
%   waterlines         - [cell array] Ordered waterline structures for all bodies, with coordinates in [m].
%   cfg                - [struct] Validated CRESTU configuration containing SI-valued physical and numerical parameters.
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

    n_body = numel(waterlines);
    outer_n = max(48, sum(cellfun(@(w)w.n_pts, waterlines)));
    theta = linspace(0, 2 * pi, outer_n + 1).';
    theta(end) = [];
    outer = cfg.fs.r_outer * [cos(theta), sin(theta)];
    points = outer;
    constraints = [(1:outer_n).', [2:outer_n, 1].'];
    offset = outer_n;
    for bodyIndex = 1:n_body
        xy = waterlines{bodyIndex}.nodes;
        nb = size(xy, 1);
        points = [points;xy]; %#ok<AGROW>
        ids = offset + (1:nb);
        constraints = [constraints;ids.', ids([2:end, 1]).']; %#ok<AGROW>
        offset = offset + nb;
    end
% Keep approximately 14 Cartesian intervals across the radius. This
% retains about 1,000 constrained triangles over two 300-panel bodies
% while the physical spacing scales with each frequency-local radius.
    spacing = max(cfg.fs.r_outer / 14, 0.75);
    axis_values = -cfg.fs.r_outer:spacing:cfg.fs.r_outer;
    [gx, gy] = meshgrid(axis_values, axis_values);
    candidate = [gx(:), gy(:)];
    keep = sum(candidate .^ 2, 2) < (0.985 * cfg.fs.r_outer)^2;
    for bodyIndex = 1:n_body
        keep = keep&~inpolygon(candidate(:, 1), candidate(:, 2), ...
            waterlines{bodyIndex}.nodes(:, 1), waterlines{bodyIndex}.nodes(:, 2));
    end
    points = [points;candidate(keep, :)];
    [points, ~, point_map] = uniquetol(points, 1e-10,'ByRows', true);
    constraints = point_map(constraints);
    constraints = constraints(constraints(:, 1) ~= constraints(:, 2), :);
    dt = delaunayTriangulation(points, constraints);
    tri = dt.ConnectivityList;
    xy = dt.Points;
    centroids = (xy(tri(:, 1), :) + xy(tri(:, 2), :) + xy(tri(:, 3), :)) / 3;
    keep = sum(centroids .^ 2, 2) <= cfg.fs.r_outer^2;
    for bodyIndex = 1:n_body
        keep = keep&~inpolygon(centroids(:, 1), centroids(:, 2), ...
            waterlines{bodyIndex}.nodes(:, 1), waterlines{bodyIndex}.nodes(:, 2));
    end
    tri = tri(keep, :);
    np = size(tri, 1);
    vertices = zeros(np, 4, 3);
    for coordinateIndex = 1:3
        vertices(:, coordinateIndex, 1:2) = xy(tri(:, coordinateIndex), :);
    end
    vertices(:, 4, :) = vertices(:, 3, :);
    centers = zeros(np, 3);
    centers(:, 1:2) = centroids(keep, :);
    normals = repmat([0, 0, 1], np, 1);
    areas = zeros(np, 1);
    e1 = zeros(np, 3);
    e2 = zeros(np, 3);
    for panelIndex = 1:np
        p1 = reshape(vertices(panelIndex, 1, :), 1, 3);
        p2 = reshape(vertices(panelIndex, 2, :), 1, 3);
        p3 = reshape(vertices(panelIndex, 3, :), 1, 3);
        areas(panelIndex) = 0.5 * norm(cross(p2 - p1, p3 - p1));
        tangent = (p2 - p1) / norm(p2 - p1);
        e1(panelIndex, :) = tangent;
        e2(panelIndex, :) = cross([0, 0, 1], tangent);
    end
    radii = sqrt(sum(centers(:, 1:2) .^ 2, 2));
    mu = zeros(np, 1);
    mask = radii > cfg.fs.r_inner;
    mu(mask) = cfg.fs.mu0 * ((radii(mask) - cfg.fs.r_inner) / (cfg.fs.r_outer - cfg.fs.r_inner)) .^ 2;
    mesh = struct('header', sprintf('Multi-body constrained FS, N=%d', n_body),'ulen', 1, ...
'panel_type', repmat(2, np, 1),'isx', 0,'isy', 0,'n_panels', np,'vertices', vertices, ...
'centers', centers,'normals', normals,'areas', areas,'e1', e1,'e2', e2, ...
'mu_damping', mu,'hydrostatics', empty_hydro());
    write_bmf(filename, mesh);
end

function hydrostatics = empty_hydro()
% EMPTY_HYDRO Return an initialized empty hydrostatic-data structure.
%
% Syntax:
%   hydrostatics = empty_hydro()
%
% Description:
%   Prepares boundary-panel geometry for the Rankine solver.
%   Global coordinates and panel-normal signs are preserved.
%
% Inputs:
%   None.
%
% Outputs:
%   hydrostatics                  - [struct] Initialized empty hydrostatic-data structure.
%
% Governing Equations / Theory:
%   Planar panel geometry, polygon moments, reflection transformations, and mesh-topology relations as applicable.
%
% References:
%   - Hess and Smith (1964); CRESTU BMF mesh-format and geometry conventions.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    hydrostatics = struct('Vx', 0,'Vy', 0,'Vz', 0,'V_mean', 0,'center_of_buoyancy', [0, 0, 0]);
end
