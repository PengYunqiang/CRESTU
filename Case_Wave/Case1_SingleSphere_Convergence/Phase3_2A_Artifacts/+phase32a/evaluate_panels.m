function result = evaluate_panels(vertices, expectedNormals, toleranceM)
% EVALUATE_PANELS Compute triangle-aware planar panel quality metrics.

    arguments
        vertices (:, 4, 3) double
        expectedNormals (:, 3) double = zeros(0, 3)
        toleranceM (1, 1) double {mustBePositive} = 1.0e-10
    end
    panelCount = size(vertices, 1);
    assert(panelCount > 0 && all(isfinite(vertices), 'all'), ...
        'CRESTU:Phase32AQualityVertices', ...
        'Panel fixtures/meshes must contain finite vertices.');
    if isempty(expectedNormals)
        expectedNormals = NaN(panelCount, 3);
    end
    assert(size(expectedNormals, 1) == panelCount, ...
        'CRESTU:Phase32AQualityNormalShape', ...
        'Expected normals must have one row per panel.');

    aspectRatio = zeros(panelCount, 1);
    minimumAngleDeg = zeros(panelCount, 1);
    maximumAngleDeg = zeros(panelCount, 1);
    skewnessDeg = zeros(panelCount, 1);
    warpDeg = zeros(panelCount, 1);
    minimumEdgeLengthM = zeros(panelCount, 1);
    areaM2 = zeros(panelCount, 1);
    windingNormal = zeros(panelCount, 3);
    orientationAlignment = NaN(panelCount, 1);
    isTriangle = false(panelCount, 1);
    unexpectedWithinPanelDuplicate = false(panelCount, 1);
    panelKeys = strings(panelCount, 1);

    for panelIndex = 1:panelCount
        panel = reshape(vertices(panelIndex, :, :), 4, 3);
        isTriangle(panelIndex) = ...
            norm(panel(4, :) - panel(3, :)) <= toleranceM;
        if isTriangle(panelIndex)
            active = panel(1:3, :);
        else
            active = panel;
        end
        quantized = int64(round(active / toleranceM));
        unexpectedWithinPanelDuplicate(panelIndex) = ...
            size(unique(quantized, 'rows'), 1) ~= size(active, 1);
        sortedVertices = sortrows(quantized);
        panelKeys(panelIndex) = string(sprintf('%d,', sortedVertices.'));

        vertexCount = size(active, 1);
        next = [2:vertexCount, 1];
        edgeVectors = active(next, :) - active;
        edgeLengths = vecnorm(edgeVectors, 2, 2);
        assert(all(edgeLengths > toleranceM), ...
            'CRESTU:Phase32AQualityZeroEdge', ...
            'A non-triangle panel contains a zero-length edge.');
        minimumEdgeLengthM(panelIndex) = min(edgeLengths);
        aspectRatio(panelIndex) = max(edgeLengths) / min(edgeLengths);

        angles = zeros(vertexCount, 1);
        for vertexIndex = 1:vertexCount
            previousIndex = mod(vertexIndex - 2, vertexCount) + 1;
            incoming = active(previousIndex, :) - active(vertexIndex, :);
            outgoing = active(next(vertexIndex), :) - active(vertexIndex, :);
            cosine = dot(incoming, outgoing) / ...
                (norm(incoming) * norm(outgoing));
            angles(vertexIndex) = acosd(max(-1, min(1, cosine)));
        end
        minimumAngleDeg(panelIndex) = min(angles);
        maximumAngleDeg(panelIndex) = max(angles);
        skewnessDeg(panelIndex) = max(abs(angles - 90));

        areaVector = cross(panel(3, :) - panel(1, :), ...
            panel(4, :) - panel(2, :));
        areaMagnitude = norm(areaVector);
        assert(areaMagnitude > 0, 'CRESTU:Phase32AQualityZeroArea', ...
            'A panel has zero area.');
        areaM2(panelIndex) = 0.5 * areaMagnitude;
        windingNormal(panelIndex, :) = areaVector / areaMagnitude;
        if all(isfinite(expectedNormals(panelIndex, :)))
            expected = expectedNormals(panelIndex, :) / ...
                norm(expectedNormals(panelIndex, :));
            orientationAlignment(panelIndex) = dot( ...
                windingNormal(panelIndex, :), expected);
        end

        if vertexCount == 4
            firstNormal = cross(edgeVectors(1, :), edgeVectors(2, :));
            secondNormal = cross(edgeVectors(3, :), edgeVectors(4, :));
            cosine = dot(firstNormal, secondNormal) / ...
                (norm(firstNormal) * norm(secondNormal));
            warpDeg(panelIndex) = acosd(max(-1, min(1, cosine)));
        end
    end

    duplicatePanelCount = panelCount - numel(unique(panelKeys));
    result = struct('schemaVersion', 1, 'panelCount', panelCount, ...
        'aspectRatio', aspectRatio, 'minimumAngleDeg', minimumAngleDeg, ...
        'maximumAngleDeg', maximumAngleDeg, 'skewnessDeg', skewnessDeg, ...
        'warpDeg', warpDeg, 'minimumEdgeLengthM', minimumEdgeLengthM, ...
        'areaM2', areaM2, 'windingNormal', windingNormal, ...
        'orientationAlignment', orientationAlignment, ...
        'isTriangle', isTriangle, ...
        'unexpectedWithinPanelDuplicate', unexpectedWithinPanelDuplicate, ...
        'duplicatePanelCount', duplicatePanelCount, ...
        'medianCharacteristicHM', median(sqrt(areaM2)), ...
        'meanCharacteristicHM', mean(sqrt(areaM2)), ...
        'totalAreaM2', sum(areaM2));
end
