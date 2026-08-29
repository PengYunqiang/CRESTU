function geometry = merge_domain_geometry(domain)
% MERGE_DOMAIN_GEOMETRY Merge domain geometry for the CRESTU hydrodynamic workflow.
%
% Syntax:
%   geometry = merge_domain_geometry(domain)
%
% Description:
%   Prepares boundary-panel geometry for the Rankine solver.
%   Global coordinates and panel-normal signs are preserved.
%
% Inputs:
%   domain             - [struct] Assembled body, free-surface, seabed, and far-field boundary domain in SI units.
%
% Outputs:
%   geometry           - [struct] Merged boundary geometry with coordinates in [m] and areas in [m^2].
%
% Governing Equations / Theory:
%   Planar panel geometry, polygon moments, reflection transformations, and mesh-topology relations as applicable.
%
% References:
%   - Hess and Smith (1964); CRESTU BMF mesh-format and geometry conventions.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    if ~isfield(domain,'body_list') || isempty(domain.body_list)
        error('CRESTU:DomainBodies','Domain has no body meshes.');
    end
    meshes = domain.body_list(:);
    names = repmat({'body'}, numel(meshes), 1);
    if isfield(domain,'fs') && ~isempty(domain.fs)
        meshes{end + 1} = domain.fs;
        names{end + 1} ='free_surface';
    end
    if isfield(domain,'seabed') && ~isempty(domain.seabed)
        meshes{end + 1} = domain.seabed;
        names{end + 1} ='seabed';
    end
    if isfield(domain,'farfield') && ~isempty(domain.farfield)
        meshes{end + 1} = domain.farfield;
        names{end + 1} ='farfield';
    end
    counts = cellfun(@(m)m.n_panels, meshes);
    total = sum(counts);
    geometry.centers = zeros(total, 3);
    geometry.normals = zeros(total, 3);
    geometry.areas = zeros(total, 1);
    geometry.vertices = zeros(total, 4, 3);
    geometry.panel_type = zeros(total, 1);
    cursor = 1;
    body_end = 0;
    for k = 1:numel(meshes)
        m = meshes{k};
        panelCount = m.n_panels;
        idx = cursor:(cursor + panelCount - 1);
        validate_mesh(m, names{k});
        geometry.centers(idx, :) = reshape(m.centers, panelCount, 3);
        geometry.normals(idx, :) = reshape(m.normals, panelCount, 3);
        geometry.areas(idx) = reshape(m.areas, panelCount, 1);
        geometry.vertices(idx, :, :) = reshape(m.vertices, panelCount, 4, 3);
        if isfield(m,'panel_type')
            geometry.panel_type(idx) = reshape(m.panel_type, panelCount, 1);
        end
        if k <= numel(domain.body_list)
            body_end = idx(end);
        end
        cursor = idx(end) + 1;
    end
    geometry.body_range = 1:body_end;
    geometry.total_panels = total;
    geometry.body_panels = body_end;
end

function validate_mesh(m, label)
% VALIDATE_MESH Validate the geometry fields required for domain merging.
%
% Syntax:
%   validate_mesh(m, label)
%
% Description:
%   Prepares boundary-panel geometry for the Rankine solver.
%   Global coordinates and panel-normal signs are preserved.
%
% Inputs:
%   m                  - [struct] Boundary mesh containing panel vertices, centers, normals, and areas in SI units.
%   label              - [character vector or string scalar] Diagnostic field label.
%
% Outputs:
%   None.
%
% Governing Equations / Theory:
%   Planar panel geometry, polygon moments, reflection transformations, and mesh-topology relations as applicable.
%
% References:
%   - Hess and Smith (1964); CRESTU BMF mesh-format and geometry conventions.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    panelCount = m.n_panels;
    required = {'centers','normals','areas','vertices'};
    for requiredFieldIndex = 1:numel(required)
        if ~isfield(m, required{requiredFieldIndex})
            error('CRESTU:MeshField','%s mesh lacks %s.', label, required{requiredFieldIndex});
        end
    end
    if numel(m.centers) ~= 3 * panelCount || numel(m.normals) ~= 3 * panelCount || numel(m.areas) ~= panelCount || numel(m.vertices) ~= 12 * panelCount
        error('CRESTU:MeshShape','%s mesh arrays do not match n_panels=%d.', label, panelCount);
    end
end
