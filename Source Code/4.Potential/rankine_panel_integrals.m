function [potentialIntegral, normalDerivativeIntegral] = ...
    rankine_panel_integrals(panelVerticesInput, fieldPoint, normal)
% RANKINE_PANEL_INTEGRALS Evaluate the analytic Rankine source and source-normal panel integrals.
%
% Syntax:
%   [potentialIntegral, normalDerivativeIntegral] = ...
%       rankine_panel_integrals(panelVerticesInput, fieldPoint, normal)
%
% Description:
%   Implements linear Rankine potential-flow operations.
%   Symmetry and phase follow exp(i*omega*t).
%
% Inputs:
%   panelVerticesInput - [4 x 3] or [3 x 4] Ordered panel vertices [m].
%   fieldPoint         - [1 x 3] Field or collocation point [m].
%   normal             - [1 x 3] Unit source-panel normal, dimensionless.
%
% Outputs:
%   potentialIntegral  - [scalar] Surface integral of the Rankine kernel 1/r [m].
%   normalDerivativeIntegral - [scalar] Source-normal derivative integral [-].
%
% Governing Equations / Theory:
%   Green third identity, the Rankine kernel 1/r, linearized free-surface theory, and reflection symmetry as applicable.
%
% References:
%   - Newman, J. N. (1977), Marine Hydrodynamics; Hess and Smith (1964); project boundary-condition specification.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    panelVertices = squeeze(panelVerticesInput);
    if size(panelVertices, 1) == 3 && size(panelVertices, 2) == 4
        panelVertices = panelVertices';
    end

    panelCenter = mean(panelVertices, 1); % [m]

    panelNormal = normal / norm(normal);
    panelTangentOne = panelVertices(2, :) - panelVertices(1, :);
    panelTangentOne = panelTangentOne - ...
        dot(panelTangentOne, panelNormal) * panelNormal;
    if norm(panelTangentOne) < 1e-12
        panelTangentOne = panelVertices(3, :) - panelVertices(2, :);
        panelTangentOne = panelTangentOne - ...
            dot(panelTangentOne, panelNormal) * panelNormal;
    end
    panelTangentOne = panelTangentOne / norm(panelTangentOne);
    panelTangentTwo = cross(panelNormal, panelTangentOne);

    transformationMatrix = [panelTangentOne; panelTangentTwo; panelNormal];

    localVertices = (panelVertices - panelCenter) * transformationMatrix';
    localFieldPoint = (fieldPoint - panelCenter) * transformationMatrix';

    localX = localFieldPoint(1);
    localY = localFieldPoint(2);
    localZ = localFieldPoint(3);

    vertexX = [localVertices(:, 1); localVertices(1, 1)];
    vertexY = [localVertices(:, 2); localVertices(1, 2)];

    potentialIntegral = 0.0;
    normalDerivativeIntegral = 0.0;
    geometricTolerance = 1e-12;

    for k = 1:4
        vertexXOne = vertexX(k);
        vertexYOne = vertexY(k);
        vertexXTwo = vertexX(k + 1);
        vertexYTwo = vertexY(k + 1);

        edgeX = vertexXTwo - vertexXOne;
        edgeY = vertexYTwo - vertexYOne;
        edgeLength = sqrt(edgeX^2 + edgeY^2); % [m]
        if edgeLength < geometricTolerance
            continue;
        end

        distanceOne = sqrt((localX - vertexXOne)^2 + ...
            (localY - vertexYOne)^2 + localZ^2); % [m]
        distanceTwo = sqrt((localX - vertexXTwo)^2 + ...
            (localY - vertexYTwo)^2 + localZ^2); % [m]

        edgeCosine = edgeX / edgeLength;
        edgeSine = edgeY / edgeLength;

        signedDistance = (vertexXOne - localX) * edgeSine - ...
            (vertexYOne - localY) * edgeCosine; % [m]

        logarithmDenominator = distanceOne + distanceTwo - edgeLength;
        if logarithmDenominator > geometricTolerance
            logarithmTerm = log((distanceOne + distanceTwo + edgeLength) / ...
                logarithmDenominator);
            potentialIntegral = potentialIntegral + signedDistance * logarithmTerm;
        end

    end

% Robust oriented solid angle.  For source-normal differentiation,
% d/dn_y(1/|X-y|) is minus the oriented angle seen from X.
    if abs(localZ) > geometricTolerance
        vertexVectorOne = panelVertices(1, :) - fieldPoint;
        vertexVectorTwo = panelVertices(2, :) - fieldPoint;
        vertexVectorThree = panelVertices(3, :) - fieldPoint;
        vertexVectorFour = panelVertices(4, :) - fieldPoint;
        solidAngle = triangle_solid_angle(vertexVectorOne, vertexVectorTwo, ...
            vertexVectorThree) + triangle_solid_angle(vertexVectorOne, ...
            vertexVectorThree, vertexVectorFour);
        normalDerivativeIntegral = -solidAngle;
    end

    potentialIntegral = potentialIntegral - localZ * normalDerivativeIntegral;

end

function solidAngle = triangle_solid_angle(vertexVectorOne, vertexVectorTwo, vertexVectorThree)
% TRIANGLE_SOLID_ANGLE Evaluate the signed solid angle subtended by a triangular facet.
%
% Syntax:
%   solidAngle = triangle_solid_angle(vertexVectorOne, vertexVectorTwo, vertexVectorThree)
%
% Description:
%   Implements linear Rankine potential-flow operations.
%   Symmetry and phase follow exp(i*omega*t).
%
% Inputs:
%   vertexVectorOne   - [1 x 3] First vertex relative to the field point [m].
%   vertexVectorTwo   - [1 x 3] Second vertex relative to the field point [m].
%   vertexVectorThree - [1 x 3] Third vertex relative to the field point [m].
%
% Outputs:
%   solidAngle         - [scalar] Signed solid angle [sr].
%
% Governing Equations / Theory:
%   Green third identity, the Rankine kernel 1/r, linearized free-surface theory, and reflection symmetry as applicable.
%
% References:
%   - Newman, J. N. (1977), Marine Hydrodynamics; Hess and Smith (1964); project boundary-condition specification.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    numerator = dot(vertexVectorOne, cross(vertexVectorTwo, vertexVectorThree));
    denominator = norm(vertexVectorOne) * norm(vertexVectorTwo) * ...
        norm(vertexVectorThree) + dot(vertexVectorOne, vertexVectorTwo) * ...
        norm(vertexVectorThree) + dot(vertexVectorTwo, vertexVectorThree) * ...
        norm(vertexVectorOne) + dot(vertexVectorThree, vertexVectorOne) * ...
        norm(vertexVectorTwo);
    solidAngle = 2 * atan2(numerator, denominator);
end
