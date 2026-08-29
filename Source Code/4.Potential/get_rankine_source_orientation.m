function orientation = get_rankine_source_orientation(varargin)
% GET_RANKINE_SOURCE_ORIENTATION Define the double-layer source convention.
%
% Syntax:
%   specification = get_rankine_source_orientation()
%   orientation = get_rankine_source_orientation(stats, normals, vertices)
%
% Description:
%   The direct Rankine identity is assembled with a canonical source normal
%   directed inward to the computational fluid domain. Existing stored mesh
%   normals are retained for physical boundary data and pressure integration.
%   This helper supplies the sole, explicit conversion applied to dG/dn_y.

    componentNames = ["body", "free_surface", "bottom", "outer_boundary"];
    componentSigns = [1, -1, 1, 1];
    convention = "double-layer-source-normal-inward-to-fluid-v1";
    signature = struct('convention', convention, ...
        'componentOrder', componentNames, ...
        'storedToCanonicalSigns', componentSigns);
    orientation = struct('convention', convention, ...
        'componentNames', componentNames, ...
        'componentSigns', componentSigns, ...
        'signatureHash', sha256_hash(jsonencode(signature)), ...
        'columnSigns', zeros(0, 1), ...
        'componentPanelCounts', zeros(1, 4), ...
        'minimumWindingAlignment', NaN(1, 4));

    if nargin == 0
        return;
    end
    if nargin ~= 3
        error('CRESTU:SourceOrientationArguments', ...
            'Expected no inputs or stats, normals, and vertices.');
    end

    stats = varargin{1};
    normals = varargin{2};
    vertices = varargin{3};
    requiredCounts = {'total_body_panels', 'fs_panels', 'seabed_panels'};
    for fieldIndex = 1:numel(requiredCounts)
        if ~isfield(stats, requiredCounts{fieldIndex})
            error('CRESTU:SourceOrientationCounts', ...
                'Panel statistics lack %s.', requiredCounts{fieldIndex});
        end
    end
    panelCount = size(normals, 1);
    normals = reshape(normals, panelCount, 3);
    vertices = reshape(vertices, panelCount, 4, 3);
    componentPanelCounts = [stats.total_body_panels, stats.fs_panels, ...
        stats.seabed_panels, panelCount - stats.total_body_panels - ...
        stats.fs_panels - stats.seabed_panels];
    if componentPanelCounts(4) < 0 || ...
            (isfield(stats, 'farfield_panels') && ...
            stats.farfield_panels ~= componentPanelCounts(4))
        error('CRESTU:SourceOrientationCounts', ...
            'Boundary-component counts do not match the geometry.');
    end

    columnSigns = zeros(panelCount, 1);
    minimumWindingAlignment = NaN(1, 4);
    firstPanel = 1;
    for componentIndex = 1:4
        count = componentPanelCounts(componentIndex);
        if count == 0
            continue;
        end
        panelRows = firstPanel:(firstPanel + count - 1);
        columnSigns(panelRows) = componentSigns(componentIndex);
        alignment = zeros(count, 1);
        for localIndex = 1:count
            panelIndex = panelRows(localIndex);
            panel = reshape(vertices(panelIndex, :, :), 4, 3);
            areaVector = cross(panel(3, :) - panel(1, :), ...
                panel(4, :) - panel(2, :));
            if norm(areaVector) <= 1e-12
                areaVector = cross(panel(2, :) - panel(1, :), ...
                    panel(4, :) - panel(1, :));
            end
            if norm(areaVector) <= 1e-12 || norm(normals(panelIndex, :)) <= 1e-12
                error('CRESTU:SourceOrientationDegeneratePanel', ...
                    'Cannot establish winding alignment for panel %d.', panelIndex);
            end
            windingNormal = areaVector / norm(areaVector);
            storedNormal = normals(panelIndex, :) / norm(normals(panelIndex, :));
            alignment(localIndex) = dot(windingNormal, storedNormal);
        end
        minimumWindingAlignment(componentIndex) = min(alignment);
        if minimumWindingAlignment(componentIndex) < 1 - 1e-10
            error('CRESTU:SourceOrientationWinding', ...
                ['%s source panels have stored-normal/winding alignment %.16g; ', ...
                'the double-layer sign conversion would be ambiguous.'], ...
                componentNames(componentIndex), ...
                minimumWindingAlignment(componentIndex));
        end
        firstPanel = panelRows(end) + 1;
    end

    orientation.columnSigns = columnSigns;
    orientation.componentPanelCounts = componentPanelCounts;
    orientation.minimumWindingAlignment = minimumWindingAlignment;
end
