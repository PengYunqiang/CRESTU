function full_waterline = complete_waterline_by_symmetry(waterline, isx, isy)
% COMPLETE_WATERLINE_BY_SYMMETRY Complete waterline by symmetry for the CRESTU hydrodynamic workflow.
%
% Syntax:
%   full_waterline = complete_waterline_by_symmetry(waterline, isx, isy)
%
% Description:
%   Prepares boundary-panel geometry for the Rankine solver.
%   Global coordinates and panel-normal signs are preserved.
%
% Inputs:
%   waterline          - [struct] Ordered waterline nodes and segment metadata, with coordinates in [m].
%   isx                - [logical scalar] Symmetry-reduction flag for the x = 0 plane.
%   isy                - [logical scalar] Symmetry-reduction flag for the y = 0 plane.
%
% Outputs:
%   full_waterline     - [struct] Full-domain waterline reconstructed from symmetry.
%
% Governing Equations / Theory:
%   Planar panel geometry, polygon moments, reflection transformations, and mesh-topology relations as applicable.
%
% References:
%   - Hess and Smith (1964); CRESTU BMF mesh-format and geometry conventions.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    points = waterline.nodes;
    if isx
        points = [points;[-points(:, 1), points(:, 2)]];
    end
    if isy
        points = [points;[points(:, 1), -points(:, 2)]];
    end
    points = uniquetol(points, 1e-9,'ByRows', true);
    center = mean(points, 1);
    angles = atan2(points(:, 2) - center(2), points(:, 1) - center(1));
    [~, order] = sort(angles);
    points = points(order, :);
    full_waterline = struct('nodes', points,'theta', atan2(points(:, 2), points(:, 1)), ...
'r', sqrt(sum(points .^ 2, 2)),'n_pts', size(points, 1),'is_closed', true);
end
