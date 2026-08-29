function generalizedNormals = compute_generalized_normals(centers, normals, body_list)
% COMPUTE_GENERALIZED_NORMALS Construct translational and rotational generalized body-normal components.
%
% Syntax:
%   generalizedNormals = compute_generalized_normals(centers, normals, body_list)
%
% Description:
%   Implements linear Rankine potential-flow operations.
%   Symmetry and phase follow exp(i*omega*t).
%
% Inputs:
%   centers            - [N x 3] Panel collocation points in global coordinates, [m].
%   normals            - [N x 3] Unit panel normals, dimensionless.
%   body_list          - [cell array] Ordered body meshes and hydrostatic metadata in SI units.
%
% Outputs:
%   generalizedNormals - [N x 6Nb] Generalized panel normals, with rotational columns in [m].
%
% Governing Equations / Theory:
%   Green third identity, the Rankine kernel 1/r, linearized free-surface theory, and reflection symmetry as applicable.
%
% References:
%   - Newman, J. N. (1977), Marine Hydrodynamics; Hess and Smith (1964); project boundary-condition specification.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    centers = normalize_xyz(centers,'centers');
    normals = normalize_xyz(normals,'normals');
    if size(centers, 1) ~= size(normals, 1)
        error('CRESTU:GeometrySize','Centers and normals must have the same row count.');
    end
    if isstruct(body_list)
        body_list = num2cell(body_list);
    end
    if ~iscell(body_list) || isempty(body_list)
        error('CRESTU:BodyList','body_list must be a nonempty cell/struct array.');
    end
    bodyCount = numel(body_list);
    panelCounts = zeros(bodyCount, 1);

    for bodyIndex = 1:bodyCount
        if ~isfield(body_list{bodyIndex}, 'n_panels') || ...
                ~isfield(body_list{bodyIndex}, 'cg')
            error('CRESTU:BodyMetadata', ...
                'Body %d requires n_panels and cg fields.', bodyIndex);
        end

        panelCounts(bodyIndex) = body_list{bodyIndex}.n_panels;
    end

    if sum(panelCounts) ~= size(centers, 1)
        error('CRESTU:BodyPanelCount', ...
            'Body panel counts total %d, geometry has %d rows.', ...
            sum(panelCounts), size(centers, 1));
    end

    %% Stage 2: Assemble translational and rotational normal components

    % A five-module model has 30 global DOFs. Module m uses
    % 6*(m-1)+(1:6). Local DOF 1-6: Surge, Sway, Heave, Roll, Pitch, Yaw.
    generalizedNormals = zeros(size(centers, 1), 6 * bodyCount);
    firstPanelIndex = 1;

    for bodyIndex = 1:bodyCount
        panelRows = firstPanelIndex:(firstPanelIndex + panelCounts(bodyIndex) - 1);
        globalDofColumns = 6 * (bodyIndex - 1) + (1:6);
        bodyNormals = normals(panelRows, :);
        centerOfGravity = reshape(body_list{bodyIndex}.cg, 1, 3); % [m]
        generalizedNormals(panelRows, globalDofColumns(1:3)) = bodyNormals;
        generalizedNormals(panelRows, globalDofColumns(4:6)) = ...
            cross(centers(panelRows, :) - centerOfGravity, bodyNormals, 2); % [m]
        firstPanelIndex = panelRows(end) + 1;
    end
end

function xyz = normalize_xyz(xyz, label)
% NORMALIZE_XYZ Normalize Cartesian coordinate storage to an N-by-3 array.
%
% Syntax:
%   xyz = normalize_xyz(xyz, label)
%
% Description:
%   Implements linear Rankine potential-flow operations.
%   Symmetry and phase follow exp(i*omega*t).
%
% Inputs:
%   xyz                - [N x 3] Cartesian coordinate array, [m].
%   label              - [character vector or string scalar] Diagnostic field label.
%
% Outputs:
%   xyz                - [N x 3] Cartesian coordinate array, [m].
%
% Governing Equations / Theory:
%   Green third identity, the Rankine kernel 1/r, linearized free-surface theory, and reflection symmetry as applicable.
%
% References:
%   - Newman, J. N. (1977), Marine Hydrodynamics; Hess and Smith (1964); project boundary-condition specification.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    if isempty(xyz)
        xyz = zeros(0, 3);
    end
    if isvector(xyz) && numel(xyz) == 3
        xyz = reshape(xyz, 1, 3);
    end
    if ~isnumeric(xyz) || ~ismatrix(xyz) || size(xyz, 2) ~= 3 || any(~isfinite(xyz(:)))
        error('CRESTU:GeometryShape','%s must be a finite N-by-3 numeric array.', label);
    end
end
