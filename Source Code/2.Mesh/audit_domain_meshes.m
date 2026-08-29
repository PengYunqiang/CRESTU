function audit = audit_domain_meshes(domain)
% AUDIT_DOMAIN_MESHES Fingerprint every boundary passed to BEM assembly.
%
% Syntax:
%   audit = audit_domain_meshes(domain)
%
% Inputs:
%   domain - Assembled body/free-surface/bottom/outer boundary domain.
%
% Outputs:
%   audit  - Component counts, SHA-256 hashes, and merged-geometry hash.

    %% 阶段 1: 验证求解域并整理物面

    assert(isstruct(domain) && isfield(domain, 'body_list') && ...
        ~isempty(domain.body_list), 'CRESTU:AuditDomainBody', ...
        'Mesh audit requires at least one body mesh.');
    assert(isfield(domain, 'geometry'), 'CRESTU:AuditDomainGeometry', ...
        'Mesh audit requires merged BEM geometry.');
    combinedBodyMesh = combine_meshes(domain.body_list);

    %% 阶段 2: 对各类边界和合并几何计算指纹

    bodyFingerprint = mesh_fingerprint(combinedBodyMesh);
    freeSurfaceFingerprint = mesh_fingerprint(optional_mesh(domain, 'fs'));
    bottomFingerprint = mesh_fingerprint(optional_mesh(domain, 'seabed'));
    outerFingerprint = mesh_fingerprint(optional_mesh(domain, 'farfield'));
    geometryMesh = geometry_as_mesh(domain.geometry);
    geometryFingerprint = mesh_fingerprint(geometryMesh);
    bodyWaterlinePanelCount = count_waterline_panels(domain, ...
        'body_waterlines', 'waterline');
    outerWaterlinePanelCount = count_waterline_panels(domain, ...
        'outer_waterlines', 'waterline');

    %% 阶段 3: 构造可保存、可比较的审计记录

    bodyFreeSurfaceHash = sha256_hash([bodyFingerprint.hash, ...
        freeSurfaceFingerprint.hash]);
    audit = struct('body', bodyFingerprint, ...
        'freeSurface', freeSurfaceFingerprint, ...
        'bottom', bottomFingerprint, ...
        'outerBoundary', outerFingerprint, ...
        'geometry', geometryFingerprint, ...
        'bodyWaterlinePanelCount', bodyWaterlinePanelCount, ...
        'outerWaterlinePanelCount', outerWaterlinePanelCount, ...
        'bodyFreeSurfaceHash', bodyFreeSurfaceHash, ...
        'bemUnknownCount', domain.geometry.total_panels, ...
        'sourcePanelCount', domain.geometry.total_panels, ...
        'collocationPointCount', size(domain.geometry.centers, 1), ...
        'sourcePointHash', geometryFingerprint.hash, ...
        'collocationPointHash', hash_numeric_array(domain.geometry.centers));
end

function mesh = combine_meshes(meshList)
% COMBINE_MESHES Concatenate body meshes without altering their coordinates.

    fields = {'vertices', 'centers', 'normals', 'areas', 'panel_type'};
    mesh = struct();

    for fieldIndex = 1:numel(fields)
        fieldName = fields{fieldIndex};
        values = cellfun(@(item) item.(fieldName), meshList, ...
            'UniformOutput', false);
        mesh.(fieldName) = cat(1, values{:});
    end

    mesh.n_panels = size(mesh.centers, 1);
end

function mesh = optional_mesh(domain, fieldName)
% OPTIONAL_MESH Return an optional boundary or an empty value.

    mesh = [];

    if isfield(domain, fieldName)
        mesh = domain.(fieldName);
    end
end

function mesh = geometry_as_mesh(geometry)
% GEOMETRY_AS_MESH Adapt merged geometry to the mesh fingerprint contract.

    mesh = struct('n_panels', geometry.total_panels, ...
        'vertices', geometry.vertices, 'centers', geometry.centers, ...
        'normals', geometry.normals, 'areas', geometry.areas, ...
        'panel_type', geometry.panel_type);
end

function panelCount = count_waterline_panels(domain, primaryField, fallbackField)
% COUNT_WATERLINE_PANELS Count actual closed/open waterline segments.

    waterlines = [];

    if isfield(domain, primaryField)
        waterlines = domain.(primaryField);
    elseif isfield(domain, fallbackField)
        waterlines = domain.(fallbackField);
    end

    if isempty(waterlines)
        panelCount = 0;
        return
    end

    if ~iscell(waterlines)
        waterlines = {waterlines};
    end

    panelCount = 0;

    for waterlineIndex = 1:numel(waterlines)
        waterline = waterlines{waterlineIndex};

        if waterline.is_closed
            panelCount = panelCount + waterline.n_pts;
        else
            panelCount = panelCount + max(0, waterline.n_pts - 1);
        end
    end
end

function hashText = hash_numeric_array(values)
% HASH_NUMERIC_ARRAY Hash dimensions and exact double-precision bytes.

    sizeBytes = reshape(typecast(uint64(size(values)), 'uint8'), 1, []);
    valueBytes = reshape(typecast(double(values(:)), 'uint8'), 1, []);
    hashText = sha256_hash([sizeBytes, valueBytes]);
end
