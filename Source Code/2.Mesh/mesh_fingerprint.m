function fingerprint = mesh_fingerprint(mesh)
% MESH_FINGERPRINT Count and hash the exact arrays used by the BEM solver.
%
% Syntax:
%   fingerprint = mesh_fingerprint(mesh)
%
% Inputs:
%   mesh        - Boundary mesh containing vertices, centers, normals, and
%                 panel areas in SI units.
%
% Outputs:
%   fingerprint - Panel/vertex counts, area, characteristic size, and
%                 deterministic SHA-256 mesh hash.

    %% 阶段 1: 处理空边界并验证数组

    if isempty(mesh)
        fingerprint = empty_fingerprint();
        return
    end

    requiredFields = {'n_panels', 'vertices', 'centers', 'normals', 'areas'};

    for fieldIndex = 1:numel(requiredFields)
        assert(isfield(mesh, requiredFields{fieldIndex}), ...
            'CRESTU:MeshFingerprintField', ...
            'Mesh selected for hashing lacks field %s.', ...
            requiredFields{fieldIndex});
    end

    panelCount = mesh.n_panels;
    vertices = reshape(double(mesh.vertices), panelCount, 4, 3); % [m]
    centers = reshape(double(mesh.centers), panelCount, 3); % [m]
    normals = reshape(double(mesh.normals), panelCount, 3); % [-]
    areas = reshape(double(mesh.areas), panelCount, 1); % [m^2]
    panelTypes = zeros(panelCount, 1);

    if isfield(mesh, 'panel_type')
        panelTypes = reshape(double(mesh.panel_type), panelCount, 1);
    end

    assert(all(isfinite(vertices(:))) && all(isfinite(centers(:))) && ...
        all(isfinite(normals(:))) && all(isfinite(areas(:))), ...
        'CRESTU:MeshFingerprintFinite', ...
        'Mesh selected for hashing contains NaN or Inf.');

    %% 阶段 2: 统计唯一顶点与特征面元尺度

    vertexRows = reshape(vertices, [], 3); % [m]
    uniqueVertices = unique(vertexRows, 'rows');
    surfaceArea = sum(areas); % [m^2]
    characteristicSize = sqrt(surfaceArea / panelCount); % [m]

    %% 阶段 3: 对实际参与求解的数组计算哈希

    canonicalBytes = [encode_size(size(vertices)), encode_double(vertices), ...
        encode_size(size(centers)), encode_double(centers), ...
        encode_size(size(normals)), encode_double(normals), ...
        encode_size(size(areas)), encode_double(areas), ...
        encode_size(size(panelTypes)), encode_double(panelTypes)];
    meshHash = sha256_hash(canonicalBytes);
    fingerprint = struct('panelCount', panelCount, ...
        'vertexCount', size(uniqueVertices, 1), ...
        'surfaceArea', surfaceArea, ...
        'characteristicSize', characteristicSize, ...
        'hash', meshHash);
end

function encodedBytes = encode_size(arraySize)
% ENCODE_SIZE Encode array dimensions without losing integer precision.

    encodedBytes = reshape(typecast(uint64(arraySize), 'uint8'), 1, []);
end

function encodedBytes = encode_double(values)
% ENCODE_DOUBLE Encode double-precision values in their exact binary form.

    encodedBytes = reshape(typecast(double(values(:)), 'uint8'), 1, []);
end

function fingerprint = empty_fingerprint()
% EMPTY_FINGERPRINT Return a stable representation of an absent boundary.

    fingerprint = struct('panelCount', 0, 'vertexCount', 0, ...
        'surfaceArea', 0.0, 'characteristicSize', NaN, ...
        'hash', sha256_hash('CRESTU_EMPTY_MESH_V1'));
end
