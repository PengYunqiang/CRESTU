function full_mesh = expand_mesh_by_symmetry(mesh, isx, isy)
% EXPAND_MESH_BY_SYMMETRY Expand mesh by symmetry for the CRESTU hydrodynamic workflow.
%
% Syntax:
%   full_mesh = expand_mesh_by_symmetry(mesh, isx, isy)
%
% Description:
%   Prepares boundary-panel geometry for the Rankine solver.
%   Global coordinates and panel-normal signs are preserved.
%
% Inputs:
%   mesh               - [struct] Boundary-panel mesh with Cartesian geometry in SI units.
%   isx                - [logical scalar] Symmetry-reduction flag for the x = 0 plane.
%   isy                - [logical scalar] Symmetry-reduction flag for the y = 0 plane.
%
% Outputs:
%   full_mesh          - [struct] Full-domain mesh reconstructed from reflection symmetry.
%
% Governing Equations / Theory:
%   Planar panel geometry, polygon moments, reflection transformations, and mesh-topology relations as applicable.
%
% References:
%   - Hess and Smith (1964); CRESTU BMF mesh-format and geometry conventions.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    full_mesh = mesh;
    if isx
        full_mesh = append_reflection(full_mesh, 1);
    end
    if isy
        full_mesh = append_reflection(full_mesh, 2);
    end
    full_mesh.isx = 0;
    full_mesh.isy = 0;
    full_mesh.header = sprintf('%s | geometric full expansion', mesh.header);
end

function out = append_reflection(in, axis_index)
% APPEND_REFLECTION Append a reflected copy of a symmetry-reduced mesh.
%
% Syntax:
%   out = append_reflection(in, axis_index)
%
% Description:
%   Prepares boundary-panel geometry for the Rankine solver.
%   Global coordinates and panel-normal signs are preserved.
%
% Inputs:
%   in                 - [character vector, string scalar, or numeric value] Value to normalize.
%   axis_index         - [scalar] Cartesian reflection-axis index, dimensionless.
%
% Outputs:
%   out                - [value] Normalized or selected result with type and physical units preserved from the input.
%
% Governing Equations / Theory:
%   Planar panel geometry, polygon moments, reflection transformations, and mesh-topology relations as applicable.
%
% References:
%   - Hess and Smith (1964); CRESTU BMF mesh-format and geometry conventions.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    reflected = in;
    reflected.vertices(:, :, axis_index) = -reflected.vertices(:, :, axis_index);
    reflected.centers(:, axis_index) = -reflected.centers(:, axis_index);
    reflected.normals(:, axis_index) = -reflected.normals(:, axis_index);
    reflected.e1(:, axis_index) = -reflected.e1(:, axis_index);
    reflected.e2(:, axis_index) = -reflected.e2(:, axis_index);
% Reflection reverses vertex orientation. Swap 2 and 4 to preserve the
% original normal represented by the reflected normal array.
    reflected.vertices(:, [2, 4], :) = reflected.vertices(:, [4, 2], :);
    out = in;
    fields = {'vertices','centers','normals','areas','e1','e2','panel_type','mu_damping'};
    for k = 1:numel(fields)
        name = fields{k};
        if isfield(in, name) && ~isempty(in.(name))
            out.(name) = cat(1, in.(name), reflected.(name));
        end
    end
    out.n_panels = size(out.centers, 1);
end
