function [fieldPotential, diagnostics] = evaluate_rankine_off_surface( ...
        fieldPoints, sourceNormals, sourceVertices, stats, ...
        boundaryPotential, canonicalNormalDerivative, cfg, parity)
% EVALUATE_RANKINE_OFF_SURFACE Reconstruct Phi in the fluid from boundary data.
%
% For canonical source normals n_s directed inward to the fluid,
%   Phi(x) = integral[dG/dn_s * Phi - G * q_s] / (4*pi),
% where q_s = dPhi/dn_s. Stored mesh normals are converted only for the
% double-layer derivative by get_rankine_source_orientation. This routine
% never changes production boundary data or the solved boundary potential.

    if nargin < 8 || isempty(parity)
        parity = [1, 1];
    end
    fieldPoints = reshape(fieldPoints, [], 3);
    panelCount = size(sourceNormals, 1);
    sourceNormals = reshape(sourceNormals, panelCount, 3);
    sourceVertices = reshape(sourceVertices, panelCount, 4, 3);
    assert(size(boundaryPotential, 1) == panelCount && ...
        isequal(size(boundaryPotential), size(canonicalNormalDerivative)), ...
        'CRESTU:OffSurfaceBoundaryShape', ...
        'Boundary potential and canonical derivative shapes must match.');
    assert(all(isfinite(fieldPoints(:))) && ...
        all(isfinite(sourceNormals(:))) && ...
        all(isfinite(sourceVertices(:))), ...
        'CRESTU:OffSurfaceFinite', ...
        'Off-surface reconstruction inputs must be finite.');
    parity = reshape(parity, 1, 2);
    assert(all(abs(parity) == 1));
    orientation = get_rankine_source_orientation( ...
        stats, sourceNormals, sourceVertices);
    [imageVertices, imageNormals, imageWeights, imageLabels] = ...
        prepare_images(sourceVertices, sourceNormals, cfg, parity);
    fieldPotential = complex(zeros(size(fieldPoints, 1), ...
        size(boundaryPotential, 2)));
    inv4pi = 1 / (4 * pi);
    for fieldIndex = 1:size(fieldPoints, 1)
        point = fieldPoints(fieldIndex, :);
        value = complex(zeros(1, size(boundaryPotential, 2)));
        for sourceIndex = 1:panelCount
            singleLayer = 0;
            canonicalDoubleLayer = 0;
            for imageIndex = 1:numel(imageWeights)
                panel = reshape(imageVertices(sourceIndex, :, :, imageIndex), ...
                    4, 3);
                storedNormal = reshape(imageNormals( ...
                    sourceIndex, :, imageIndex), 1, 3);
                [greenIntegral, storedDerivativeIntegral] = ...
                    rankine_panel_integrals(panel, point, storedNormal);
                weight = imageWeights(imageIndex);
                singleLayer = singleLayer + weight * greenIntegral;
                canonicalDoubleLayer = canonicalDoubleLayer + weight * ...
                    orientation.columnSigns(sourceIndex) * ...
                    storedDerivativeIntegral;
            end
            value = value + inv4pi * (canonicalDoubleLayer * ...
                boundaryPotential(sourceIndex, :) - singleLayer * ...
                canonicalNormalDerivative(sourceIndex, :));
        end
        fieldPotential(fieldIndex, :) = value;
    end
    diagnostics = struct('schemaVersion', 1, ...
        'timeConvention', 'exp(+i*omega*t)', ...
        'representation', ...
        'Phi=int[(dG/dn_s)Phi-G*q_s]dS/(4*pi)', ...
        'canonicalNormalConvention', char(orientation.convention), ...
        'orientationSHA256', orientation.signatureHash, ...
        'fieldPointCount', size(fieldPoints, 1), ...
        'sourcePanelCount', panelCount, ...
        'columnCount', size(boundaryPotential, 2), ...
        'parity', parity, 'imageLabels', {imageLabels});
end

function [allVertices, allNormals, weights, labels] = ...
        prepare_images(vertices, normals, cfg, parity)
% PREPARE_IMAGES Match the production assembler's reflection convention.
    flags = [0, 0];
    if cfg.isx, flags = [flags; 1, 0]; end
    if cfg.isy, flags = [flags; 0, 1]; end
    if cfg.isx && cfg.isy, flags = [flags; 1, 1]; end
    panelCount = size(vertices, 1);
    imageCount = size(flags, 1);
    allVertices = zeros(panelCount, 4, 3, imageCount);
    allNormals = zeros(panelCount, 3, imageCount);
    weights = ones(imageCount, 1);
    labels = cell(imageCount, 1);
    for imageIndex = 1:imageCount
        imageVertices = vertices;
        imageNormal = normals;
        if flags(imageIndex, 1)
            imageVertices(:, :, 1) = -imageVertices(:, :, 1);
            imageNormal(:, 1) = -imageNormal(:, 1);
        end
        if flags(imageIndex, 2)
            imageVertices(:, :, 2) = -imageVertices(:, :, 2);
            imageNormal(:, 2) = -imageNormal(:, 2);
        end
        if mod(sum(flags(imageIndex, :)), 2) == 1
            imageVertices(:, [2, 4], :) = ...
                imageVertices(:, [4, 2], :);
        end
        allVertices(:, :, :, imageIndex) = imageVertices;
        allNormals(:, :, imageIndex) = imageNormal;
        weights(imageIndex) = parity(1)^flags(imageIndex, 1) * ...
            parity(2)^flags(imageIndex, 2);
        labels{imageIndex} = sprintf('reflect_x=%d,reflect_y=%d', ...
            flags(imageIndex, 1), flags(imageIndex, 2));
    end
end
