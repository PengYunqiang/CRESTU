function result = solve_uniform_flow(mesh_body, U_inf)
% SOLVE_UNIFORM_FLOW Solve the Hess-Smith source-panel system for steady uniform flow.
%
% Syntax:
%   result = solve_uniform_flow(mesh_body, U_inf)
%
% Description:
%   Solve the Hess-Smith source-panel system for steady uniform flow.
%   The implementation preserves the CRESTU solver and data-field contracts.
%
% Inputs:
%   mesh_body - Body-panel mesh structure with coordinates in meters [m].
%   U_inf - Free-stream velocity vector, [1 x 3] [m/s].
%
% Outputs:
%   result - Computed hydrodynamic result structure [-].
%
% Governing Equations / Theory:
%   Uses the linear potential-flow, source-panel, or post-processing
%   formulation stated in the CRESTU theory and technical manual.
%
% References:
%   - CRESTU theory and technical manual.
%   - Standard boundary-element and marine hydrodynamics references.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)
    arguments
        mesh_body (1, 1) struct
        U_inf (1, 3) double {mustBeFinite}
    end

    %% Stage 1: Validate inputs and geometry

    requiredFields = {'n_panels','centers','normals','vertices','e1','e2'};

    for fieldIndex = 1:numel(requiredFields)
        assert(isfield(mesh_body, requiredFields{fieldIndex}), ...
'CRESTU:UniformFlowMeshField', ...
'Body mesh is missing field "%s".', requiredFields{fieldIndex});
    end

    panelCount = mesh_body.n_panels;
    assert(panelCount > 0,'CRESTU:UniformFlowEmptyMesh', ...
'The body mesh must contain at least one panel.');
    assert(isequal(size(mesh_body.centers), [panelCount, 3]), ...
'CRESTU:UniformFlowCenterShape', ...
'Body centers must have size [N_panels x 3].');
    assert(isequal(size(mesh_body.normals), [panelCount, 3]), ...
'CRESTU:UniformFlowNormalShape', ...
'Body normals must have size [N_panels x 3].');
    assert(size(mesh_body.vertices, 1) == panelCount && ...
        size(mesh_body.vertices, 2) == 4 && size(mesh_body.vertices, 3) == 3, ...
'CRESTU:UniformFlowVertexShape', ...
'Body vertices must have size [N_panels x 4 x 3].');
    assert(dot(U_inf, U_inf) > 0,'CRESTU:UniformFlowZeroSpeed', ...
'Free-stream speed must be positive.');

    panelCenters = mesh_body.centers; % [m]
    panelNormals = mesh_body.normals; % [-]
    panelVertices = mesh_body.vertices; % [m]
    panelTangentOne = mesh_body.e1; % [-]
    panelTangentTwo = mesh_body.e2; % [-]

    %% Stage 2: Assemble the Hess-Smith influence matrices

    fprintf('[INFO] Assemble the BEM influence matrix | panels = %d\n', panelCount);
    assemblyTimer = tic;

    A = zeros(panelCount, panelCount); % Normal-velocity influence matrix [-]
    boundaryRhs = -panelNormals * U_inf(:); % Boundary normal velocity [m/s]
    velocityInfluenceX = zeros(panelCount, panelCount); % [-]
    velocityInfluenceY = zeros(panelCount, panelCount); % [-]
    velocityInfluenceZ = zeros(panelCount, panelCount); % [-]

    parfor i = 1:panelCount
        fieldPoint = panelCenters(i, :); % [m]
        fieldNormal = panelNormals(i, :); % [-]
        normalInfluenceRow = zeros(1, panelCount);
        velocityInfluenceXRow = zeros(1, panelCount);
        velocityInfluenceYRow = zeros(1, panelCount);
        velocityInfluenceZRow = zeros(1, panelCount);

        normalInfluenceRow(i) = 0.5; % Self-influence jump term.

        for j = 1:panelCount
            if i ~= j
                sourcePanelVertices = squeeze(panelVertices(j, :, :)); % [m]
                [inducedVelocityX, inducedVelocityY, inducedVelocityZ] = ...
                    hess_smith_panel_velocity(sourcePanelVertices, fieldPoint, ...
                    panelCenters(j, :), panelTangentOne(j, :), ...
                    panelTangentTwo(j, :), panelNormals(j, :));

                velocityInfluenceXRow(j) = inducedVelocityX;
                velocityInfluenceYRow(j) = inducedVelocityY;
                velocityInfluenceZRow(j) = inducedVelocityZ;
                normalInfluenceRow(j) = inducedVelocityX * fieldNormal(1) + ...
                    inducedVelocityY * fieldNormal(2) + ...
                    inducedVelocityZ * fieldNormal(3);
            end
        end

        A(i, :) = normalInfluenceRow;
        velocityInfluenceX(i, :) = velocityInfluenceXRow;
        velocityInfluenceY(i, :) = velocityInfluenceYRow;
        velocityInfluenceZ(i, :) = velocityInfluenceZRow;
    end

    fprintf('[OK] Influence matrix assembled | elapsed = %.2f s\n', toc(assemblyTimer));

    %% Stage 3: Solve source strengths and recover surface velocity

    sourceStrength = A \ boundaryRhs; % [m/s]
    inducedVelocity = [velocityInfluenceX * sourceStrength, ...
        velocityInfluenceY * sourceStrength, ...
        velocityInfluenceZ * sourceStrength]; % [m / s]
    totalVelocity = repmat(U_inf, panelCount, 1) + inducedVelocity + ...
        0.5 * (sourceStrength .* panelNormals); % [m/s]
    freeStreamSpeedSquared = dot(U_inf, U_inf); % [m^2/s^2]
    pressureCoefficient = 1.0 - ...
        sum(totalVelocity .^ 2, 2) / freeStreamSpeedSquared; % [-]

    %% Stage 4: Preserve the public result-field contract

    result = struct();
    result.sigma = sourceStrength;
    result.V_total = totalVelocity;
    result.Cp = pressureCoefficient;

    fprintf('[OK] Uniform-flow solution completed | panels = %d\n', panelCount);
end
