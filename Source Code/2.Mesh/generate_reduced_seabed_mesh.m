function mesh = generate_reduced_seabed_mesh(cfg, full_circumference_nodes)
% GENERATE_REDUCED_SEABED_MESH Generate reduced seabed mesh for the CRESTU hydrodynamic workflow.
%
% Syntax:
%   mesh = generate_reduced_seabed_mesh(cfg, full_circumference_nodes)
%
% Description:
%   Prepares boundary-panel geometry for the Rankine solver.
%   Global coordinates and panel-normal signs are preserved.
%
% Inputs:
%   cfg                - [struct] Validated CRESTU configuration containing SI-valued physical and numerical parameters.
%   full_circumference_nodes - [N x 2] Full outer-boundary nodes in counterclockwise order, [m].
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

    if cfg.water_depth <= 0
        error('CRESTU:SeabedDepth','Finite positive depth required.');
    end
    if cfg.isx && cfg.isy
        theta0 = 0;
        theta1 = pi / 2;
        divisor = 4;
    elseif cfg.isx
        theta0 = -pi / 2;
        theta1 = pi / 2;
        divisor = 2;
    elseif cfg.isy
        theta0 = 0;
        theta1 = pi;
        divisor = 2;
    else, theta0 = 0;
    theta1 = 2 * pi;
    divisor = 1;
    end
    ntheta = max(4, round(full_circumference_nodes / divisor));
    theta = linspace(theta0, theta1, ntheta + 1);
    nr = max(4, cfg.fs.nr_near + cfg.fs.nr_sponge);
    radial = cfg.fs.r_outer * (linspace(0, 1, nr + 1) .^ 0.8);
    np = nr * ntheta;
    vertices = zeros(np, 4, 3);
    panelIndex = 0;
    z = -cfg.water_depth;
    for i = 1:nr
        for j = 1:ntheta
            panelIndex = panelIndex + 1;
            r0 = radial(i);
            r1 = radial(i + 1);
            t0 = theta(j);
            t1 = theta(j + 1);
            vertices(panelIndex, 1, :) = [r0 * cos(t0), r0 * sin(t0), z];
            vertices(panelIndex, 2, :) = [r1 * cos(t0), r1 * sin(t0), z];
            vertices(panelIndex, 3, :) = [r1 * cos(t1), r1 * sin(t1), z];
            vertices(panelIndex, 4, :) = [r0 * cos(t1), r0 * sin(t1), z];
        end
    end
    [centers, areas, e1, e2] = panel_geometry(vertices);
    mesh = struct('header', sprintf('Reduced polar seabed ISX=%d ISY=%d', cfg.isx, cfg.isy), ...
'ulen', 1,'panel_type', repmat(4, np, 1),'isx', cfg.isx,'isy', cfg.isy, ...
'n_panels', np,'vertices', vertices,'centers', centers,'normals', repmat([0, 0, 1], np, 1), ...
'areas', areas,'e1', e1,'e2', e2,'hydrostatics', empty_hydro());
end

function [centers, areas, e1, e2] = panel_geometry(vertices)
% PANEL_GEOMETRY Evaluate panel centers, areas, and local tangent bases.
%
% Syntax:
%   [centers, areas, e1, e2] = panel_geometry(vertices)
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
        panelVertices = reshape(vertices(panelIndex, :, :), 4, 3);
        centers(panelIndex, :) = mean(panelVertices, 1);
        areas(panelIndex) = 0.5 * norm(cross(...
            panelVertices(3, :) - panelVertices(1, :), ...
            panelVertices(4, :) - panelVertices(2, :)));
        tangent = panelVertices(2, :) - panelVertices(1, :);
        if norm(tangent) < 1e-12
            tangent = panelVertices(3, :) - panelVertices(1, :);
        end
        e1(panelIndex, :) = tangent / norm(tangent);
        e2(panelIndex, :) = cross([0, 0, 1], e1(panelIndex, :));
    end
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
