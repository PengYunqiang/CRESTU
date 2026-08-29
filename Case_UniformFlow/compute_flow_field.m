function [X_grid, Y_grid, Z_grid, V_field, Cp_field] = compute_flow_field(mesh_body, sigma, U_inf, domain_lim, N_pts)
% COMPUTE_FLOW_FIELD Evaluate the velocity and pressure-coefficient fields around a source-panel body.
%
% Syntax:
%   [X_grid, Y_grid, Z_grid, V_field, Cp_field] = compute_flow_field(mesh_body, sigma, U_inf, domain_lim, N_pts)
%
% Description:
%   Evaluate the velocity and pressure-coefficient fields around a source-panel body.
%   The implementation preserves the CRESTU solver and data-field contracts.
%
% Inputs:
%   mesh_body - Body-panel mesh structure with coordinates in meters [m].
%   sigma - Constant panel source strengths, [N_panels x 1] [m/s].
%   U_inf - Free-stream velocity vector, [1 x 3] [m/s].
%   domain_lim - Symmetric Cartesian half-width, positive scalar [m].
%   N_pts - Number of field points per coordinate direction, positive integer [-].
%
% Outputs:
%   X_grid - Cartesian x-coordinate grid [m].
%   Y_grid - Cartesian y-coordinate grid [m].
%   Z_grid - Cartesian z-coordinate grid [m].
%   V_field - Velocity-magnitude field [m/s].
%   Cp_field - Pressure-coefficient field [-].
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
        sigma (:, 1) double {mustBeFinite}
        U_inf (1, 3) double {mustBeFinite}
        domain_lim (1, 1) double {mustBePositive, mustBeFinite}
        N_pts (1, 1) double {mustBeInteger, mustBeGreaterThan(N_pts, 1)}
    end

    %% Stage 1: Validate inputs and cache panel geometry

    requiredFields = {'n_panels','centers','normals','vertices','e1','e2'};

    for fieldIndex = 1:numel(requiredFields)
        assert(isfield(mesh_body, requiredFields{fieldIndex}), ...
'CRESTU:FlowFieldMeshField', ...
'Body mesh is missing field "%s".', requiredFields{fieldIndex});
    end

    panelCount = mesh_body.n_panels;
    assert(numel(sigma) == panelCount,'CRESTU:FlowFieldSourceShape', ...
'Source strength must contain one value per panel.');
    assert(dot(U_inf, U_inf) > 0,'CRESTU:FlowFieldZeroSpeed', ...
'Free-stream speed must be positive.');

    panelCenters = mesh_body.centers; % [m]
    panelNormals = mesh_body.normals; % [-]
    panelVertices = mesh_body.vertices; % [m]
    panelTangentOne = mesh_body.e1; % [-]
    panelTangentTwo = mesh_body.e2; % [-]

    %% Stage 2: Build the Cartesian evaluation grid

    xCoordinates = linspace(-domain_lim, domain_lim, N_pts); % [m]
    yCoordinates = linspace(-domain_lim, domain_lim, N_pts); % [m]
    zCoordinates = linspace(-domain_lim, domain_lim, N_pts); % [m]
    [X_grid, Y_grid, Z_grid] = ndgrid(xCoordinates, yCoordinates, zCoordinates);

    fieldPointCount = numel(X_grid);
    fieldPoints = [X_grid(:), Y_grid(:), Z_grid(:)]; % [m]
    inducedVelocity = zeros(fieldPointCount, 3); % [m/s]
    freeStreamSpeedSquared = dot(U_inf, U_inf); % [m^2/s^2]

    %% Stage 3: Evaluate the panel-induced velocity field

    fprintf('[INFO] Evaluate the external flow field | points = %d\n', fieldPointCount);
    fieldTimer = tic;

    parfor i = 1:fieldPointCount
        fieldPoint = fieldPoints(i, :); % [m]
        pointVelocity = [0, 0, 0]; % [m/s]

        for j = 1:panelCount
            sourcePanelVertices = squeeze(panelVertices(j, :, :)); % [m]
            [inducedVelocityX, inducedVelocityY, inducedVelocityZ] = ...
                hess_smith_panel_velocity(sourcePanelVertices, fieldPoint, ...
                panelCenters(j, :), panelTangentOne(j, :), ...
                panelTangentTwo(j, :), panelNormals(j, :));
            pointVelocity = pointVelocity + sigma(j) * ...
                [inducedVelocityX, inducedVelocityY, inducedVelocityZ];
        end

        inducedVelocity(i, :) = pointVelocity;
    end

    fprintf('[OK] External flow field completed | elapsed = %.2f s\n', toc(fieldTimer));

    %% Stage 4: Mask the body interior and assemble output grids

    totalVelocity = repmat(U_inf, fieldPointCount, 1) + inducedVelocity; % [m/s]
    radialDistance = sqrt(sum(fieldPoints .^ 2, 2)); % [m]
    bodyRadiusEstimate = max(sqrt(sum(panelCenters .^ 2, 2))); % [m]
    insideBody = radialDistance < 0.96 * bodyRadiusEstimate;
    totalVelocity(insideBody, :) = NaN;

    pointPressureCoefficient = 1.0 - ...
        sum(totalVelocity .^ 2, 2) / freeStreamSpeedSquared; % [-]
    V_field = struct();
    V_field.u = reshape(totalVelocity(:, 1), size(X_grid));
    V_field.v = reshape(totalVelocity(:, 2), size(X_grid));
    V_field.w = reshape(totalVelocity(:, 3), size(X_grid));
    V_field.mag = reshape(sqrt(sum(totalVelocity .^ 2, 2)), size(X_grid));
    Cp_field = reshape(pointPressureCoefficient, size(X_grid));
end

% =========================================================================
% 2. Plot orthogonal field slices with streamlines and the body section
% =========================================================================
function plot_field_slices(X, Y, Z, V, Cp, mesh_body, diameter)
% PLOT_FIELD_SLICES Plot velocity and pressure fields on orthogonal planes through the body.
%
% Syntax:
%   plot_field_slices(X, Y, Z, V, Cp, mesh_body, diameter)
%
% Description:
%   Plot velocity and pressure fields on orthogonal planes through the body.
%   The implementation preserves the CRESTU solver and data-field contracts.
%
% Inputs:
%   X - Cartesian x-coordinate grid [m].
%   Y - Cartesian y-coordinate grid [m].
%   Z - Cartesian z-coordinate grid [m].
%   V - Velocity-magnitude grid [m/s].
%   Cp - Pressure-coefficient grid [-].
%   mesh_body - Body-panel mesh structure with coordinates in meters [m].
%   diameter - Body diameter [m].
%
% Outputs:
%   None.
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
%% Stage 1: Initialize inputs and dependencies

    fig = figure('Color','w','Position', [150, 100, 1100, 500],'Name','Flow Field Slices');

% (1) Y = 0 longitudinal symmetry-plane velocity field
    subplot(1, 2, 1);
    hold on;
    grid on;
    axis equal;

%% Stage 2: Run the core calculation

    mid_y = round(size(Y, 2) / 2);
    x_sub = squeeze(X(:, mid_y, :));
    z_sub = squeeze(Z(:, mid_y, :));
    mag_sub = squeeze(V.mag(:, mid_y, :));
    u_sub = squeeze(V.u(:, mid_y, :));
    w_sub = squeeze(V.w(:, mid_y, :));

    pcolor(x_sub, z_sub, mag_sub);
    shading interp;
    colormap(gca, jet);
    colorbarHandle = colorbar;
    colorbarHandle.Label.String ='|V| / U_{\infty}';

% Add streamlines.
    [sx, sz] = meshgrid(linspace(min(x_sub(:)), min(x_sub(:)), 15), linspace(min(z_sub(:)), max(z_sub(:)), 15));
    streamline(x_sub', z_sub', u_sub', w_sub', sx, sz);

% Draw the sphere section.
    th = linspace(0, 2 * pi, 100);
    fill(0.5 * diameter * cos(th), 0.5 * diameter * sin(th), [0.8 0.8 0.8],'EdgeColor','k','LineWidth', 1.2);

    title('Y = 0 longitudinal-plane velocity magnitude and streamlines','FontSize', 11,'FontWeight','bold');
    xlabel('X (m)');
    ylabel('Z (m)');
    xlim([-diameter, diameter]);
    ylim([-diameter, diameter]);

% (2) Z = 0 horizontal equatorial-plane velocity field
    subplot(1, 2, 2);
    hold on;
    grid on;
    axis equal;
    mid_z = round(size(Z, 3) / 2);
    x_sub2 = squeeze(X(:, :, mid_z));
    y_sub2 = squeeze(Y(:, :, mid_z));
    mag_sub2 = squeeze(V.mag(:, :, mid_z));
    u_sub2 = squeeze(V.u(:, :, mid_z));
    v_sub2 = squeeze(V.v(:, :, mid_z));

    pcolor(x_sub2, y_sub2, mag_sub2);
    shading interp;
    colormap(gca, jet);
    c2 = colorbar;
    c2.Label.String ='|V| / U_{\infty}';

    [sx2, sy2] = meshgrid(linspace(min(x_sub2(:)), min(x_sub2(:)), 15), linspace(min(y_sub2(:)), max(y_sub2(:)), 15));
    streamline(x_sub2', y_sub2', u_sub2', v_sub2', sx2, sy2);
    fill(0.5 * diameter * cos(th), 0.5 * diameter * sin(th), [0.8 0.8 0.8],'EdgeColor','k','LineWidth', 1.2);

    title('Z = 0 horizontal-plane velocity magnitude and streamlines','FontSize', 11,'FontWeight','bold');
    xlabel('X (m)');
    ylabel('Y (m)');
    xlim([-diameter, diameter]);
    ylim([-diameter, diameter]);
end
