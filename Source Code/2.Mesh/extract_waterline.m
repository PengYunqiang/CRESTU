function waterline = extract_waterline(mesh_body, z_tol)
% EXTRACT_WATERLINE Extract waterline for the CRESTU hydrodynamic workflow.
%
% Syntax:
%   waterline = extract_waterline(mesh_body, z_tol)
%
% Description:
%   Prepares boundary-panel geometry for the Rankine solver.
%   Global coordinates and panel-normal signs are preserved.
%
% Inputs:
%   mesh_body          - [struct] Wetted-body panel mesh with Cartesian geometry in SI units.
%   z_tol              - [scalar] Waterline elevation tolerance, [m].
%
% Outputs:
%   waterline          - [struct] Ordered closed waterline geometry with coordinates in [m].
%
% Governing Equations / Theory:
%   Planar panel geometry, polygon moments, reflection transformations, and mesh-topology relations as applicable.
%
% References:
%   - Hess and Smith (1964); CRESTU BMF mesh-format and geometry conventions.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    if nargin < 2
        z_tol = 1e-3;
    end

    verts = mesh_body.vertices;
    normals = mesh_body.normals;
    np = size(verts, 1);

    edges = [];
    for k = 1:np
        if norm(normals(k, 1:2)) < 0.2
            continue;
        end

        panelVertices = squeeze(verts(k, :, :));
        is_top = abs(panelVertices(:, 3)) <= z_tol;
        idx_top = find(is_top);

        if length(idx_top) == 2
            edges = [edges; panelVertices(idx_top(1), 1:2), panelVertices(idx_top(2), 1:2)]; %#ok<AGROW>
        end
    end

    if isempty(edges)
        error('No exterior waterline could be identified in the body mesh!');
    end

    tol_dist = 1e-3;
    n_edges = size(edges, 1);
    used = false(n_edges, 1);

    chain = [edges(1, 1:2); edges(1, 3:4)];
    used(1) = true;
% Extend both ends.  A one-direction walk truncates an open half/quarter
% waterline whenever the arbitrary first edge lies inside the chain.
    while any(~used)
        found = false;
        for i = 1:n_edges
            if used(i)
                continue;
            end
            p_start = edges(i, 1:2);
            p_end = edges(i, 3:4);
            if norm(chain(end, :) - p_start) < tol_dist
                chain = [chain;p_end];
                used(i) = true;
                found = true;
                break; %#ok<AGROW>
            elseif norm(chain(end, :) - p_end) < tol_dist
                chain = [chain;p_start];
                used(i) = true;
                found = true;
                break; %#ok<AGROW>
            elseif norm(chain(1, :) - p_end) < tol_dist
                chain = [p_start;chain];
                used(i) = true;
                found = true;
                break; %#ok<AGROW>
            elseif norm(chain(1, :) - p_start) < tol_dist
                chain = [p_end;chain];
                used(i) = true;
                found = true;
                break; %#ok<AGROW>
            end
        end
        if ~found
            break;
        end
    end

    unique_chain = chain(1, :);
    for i = 2:size(chain, 1)
        if norm(chain(i, :) - unique_chain(end, :)) > tol_dist
            unique_chain = [unique_chain; chain(i, :)]; %#ok<AGROW>
        end
    end

    is_closed = norm(unique_chain(1, :) - unique_chain(end, :)) < tol_dist || ...
                (mesh_body.isx == 0 && mesh_body.isy == 0);

    if is_closed && norm(unique_chain(1, :) - unique_chain(end, :)) < tol_dist
        unique_chain(end, :) = [];
    end

    x = unique_chain(:, 1);
    y = unique_chain(:, 2);
    area2 = sum(x .* [y(2:end); y(1)] - y .* [x(2:end); x(1)]);
    if is_closed && area2 < 0
        unique_chain = flipud(unique_chain);
    end

    waterline.nodes = unique_chain;
    waterline.theta = atan2(unique_chain(:, 2), unique_chain(:, 1));
    waterline.r = sqrt(unique_chain(:, 1) .^ 2 + unique_chain(:, 2) .^ 2);
    waterline.n_pts = size(unique_chain, 1);
    waterline.is_closed = is_closed;

    fprintf('[OK] Waterline extraction completed: with %d continuous nodes (closed: %d)\n', waterline.n_pts, is_closed);
end
