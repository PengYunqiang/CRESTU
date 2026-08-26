function nj = compute_generalized_normals(centers, normals, body_list)
% COMPUTE_GENERALIZED_NORMALS Construct translational and rotational generalized body-normal components.
%
% Syntax:
%   nj = compute_generalized_normals(centers, normals, body_list)
%
% Description:
%   The routine implements a component of the linear Rankine boundary-element formulation for incompressible, irrotational gravity-wave flow. Geometry, reflection parity, free-surface impedance, and complex phase follow the project convention exp(i*omega*t).
%
% Inputs:
%   centers            - [N x 3] Panel collocation points in global coordinates, [m].
%   normals            - [N x 3] Unit panel normals, dimensionless.
%   body_list          - [cell array] Ordered body meshes and hydrostatic metadata in SI units.
%
% Outputs:
%   nj                 - [N x 6Nb] Generalized panel normals, with rotational columns in [m].
%
% Governing Equations / Theory:
%   Green third identity, the Rankine kernel 1/r, linearized free-surface theory, and reflection symmetry as applicable.
%
% References:
%   - Newman, J. N. (1977), Marine Hydrodynamics; Hess and Smith (1964); project boundary-condition specification.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% --- 1. Validate Inputs and Initialize the Algorithm ---

    centers = normalize_xyz(centers, 'centers'); normals = normalize_xyz(normals, 'normals');
    if size(centers, 1) ~= size(normals, 1)
        error('CRESTU:GeometrySize', 'Centers and normals must have the same row count.');
    end
    if isstruct(body_list), body_list = num2cell(body_list); end
    if ~iscell(body_list) || isempty(body_list)
        error('CRESTU:BodyList', 'body_list must be a nonempty cell/struct array.');
    end
    n_bodies = numel(body_list); counts = zeros(n_bodies, 1);
    for b = 1:n_bodies
        if ~isfield(body_list{b}, 'n_panels') || ~isfield(body_list{b}, 'cg')
            error('CRESTU:BodyMetadata', 'Body %d requires n_panels and cg fields.', b);
        end
        counts(b) = body_list{b}.n_panels;
    end
    if sum(counts) ~= size(centers, 1)
        error('CRESTU:BodyPanelCount', 'Body panel counts total %d, geometry has %d rows.', ...
            sum(counts), size(centers, 1));
    end
    nj = zeros(size(centers, 1), 6 * n_bodies); first = 1;
    for b = 1:n_bodies
        idx = first:(first + counts(b) - 1); cols = (b - 1) * 6 + (1:6);
        n_body = normals(idx, :); cg = reshape(body_list{b}.cg, 1, 3);
        nj(idx, cols(1:3)) = n_body;
        nj(idx, cols(4:6)) = cross(centers(idx, :) - cg, n_body, 2);
        first = idx(end) + 1;
    end
end

function xyz = normalize_xyz(xyz, label)
% NORMALIZE_XYZ Normalize Cartesian coordinate storage to an N-by-3 array.
%
% Syntax:
%   xyz = normalize_xyz(xyz, label)
%
% Description:
%   The routine implements a component of the linear Rankine boundary-element formulation for incompressible, irrotational gravity-wave flow. Geometry, reflection parity, free-surface impedance, and complex phase follow the project convention exp(i*omega*t).
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

%% --- 1. Validate Inputs and Initialize the Algorithm ---

    if isempty(xyz), xyz = zeros(0, 3); end
    if isvector(xyz) && numel(xyz) == 3, xyz = reshape(xyz, 1, 3); end
    if ~isnumeric(xyz) || ~ismatrix(xyz) || size(xyz, 2) ~= 3 || any(~isfinite(xyz(:)))
        error('CRESTU:GeometryShape', '%s must be a finite N-by-3 numeric array.', label);
    end
end
