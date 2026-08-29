function plot_field_slices(xGrid, yGrid, zGrid, velocityField, pressureCoefficientGrid, bodyMesh, diameter)
% PLOT_FIELD_SLICES Plot velocity fields on orthogonal planes through a body.
%
% Syntax:
%   plot_field_slices(xGrid, yGrid, zGrid, velocityField, pressureCoefficientGrid, bodyMesh, diameter)
%
% Description:
%   Plots velocity magnitude and streamlines on the y = 0 and z = 0 planes.
%   The pressure grid and body mesh are checked for interface consistency.
%
% Inputs:
%   xGrid                   - Cartesian x-coordinate grid [m].
%   yGrid                   - Cartesian y-coordinate grid [m].
%   zGrid                   - Cartesian z-coordinate grid [m].
%   velocityField           - Velocity components and magnitude [m/s].
%   pressureCoefficientGrid - Pressure-coefficient grid [-].
%   bodyMesh                - Body-panel mesh structure with coordinates [m].
%   diameter                - Body diameter [m].
%
% Outputs:
%   None.
%
% Governing Equations / Theory:
%   Streamlines follow dx/ds = u and dz/ds = w or dx/ds = u and dy/ds = v.
%
% References:
%   - MATLAB graphics documentation.
%   - CRESTU theory and technical manual.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

    arguments
        xGrid double
        yGrid double
        zGrid double
        velocityField (1, 1) struct
        pressureCoefficientGrid double
        bodyMesh (1, 1) struct
        diameter (1, 1) double {mustBePositive, mustBeFinite}
    end

    %% Stage 1: Validate grid and field dimensions

    assert(isequal(size(xGrid), size(yGrid), size(zGrid), ...
        size(pressureCoefficientGrid)),'CRESTU:FlowSliceGridShape', ...
'Coordinate and pressure grids must have identical dimensions.');
    requiredVelocityFields = {'u','v','w','mag'};

    for fieldIndex = 1:numel(requiredVelocityFields)
        fieldName = requiredVelocityFields{fieldIndex};
        assert(isfield(velocityField, fieldName) && ...
            isequal(size(velocityField.(fieldName)), size(xGrid)), ...
'CRESTU:FlowSliceVelocityShape', ...
'Velocity field "%s" must match the coordinate grid.', fieldName);
    end

    assert(isfield(bodyMesh,'centers') && ~isempty(bodyMesh.centers), ...
'CRESTU:FlowSliceBodyMesh','Body mesh centers are required.');

    figure('Color','w','Position', [150, 100, 1100, 500], ...
'Name','Flow Field Slices');

    %% Stage 2: Plot the y = 0 longitudinal plane

    subplot(1, 2, 1);
    hold on;
    grid on;
    axis equal;

    middleYIndex = round(size(yGrid, 2) / 2);
    longitudinalX = squeeze(xGrid(:, middleYIndex, :)); % [m]
    longitudinalZ = squeeze(zGrid(:, middleYIndex, :)); % [m]
    longitudinalSpeed = squeeze(velocityField.mag(:, middleYIndex, :)); % [m/s]
    longitudinalVelocityX = squeeze(velocityField.u(:, middleYIndex, :)); % [m/s]
    longitudinalVelocityZ = squeeze(velocityField.w(:, middleYIndex, :)); % [m/s]

    pcolor(longitudinalX, longitudinalZ, longitudinalSpeed);
    shading interp;
    colormap(gca, jet);
    colorbarHandle = colorbar;
    colorbarHandle.Label.String ='|V| / U_{\infty}';

    [streamlineStartX, streamlineStartZ] = meshgrid(...
        min(longitudinalX(:)), ...
        linspace(min(longitudinalZ(:)), max(longitudinalZ(:)), 15));
    streamline(longitudinalX', longitudinalZ', longitudinalVelocityX', ...
        longitudinalVelocityZ', streamlineStartX, streamlineStartZ);

    sectionAngle = linspace(0, 2 * pi, 100); % [rad]
    fill(0.5 * diameter * cos(sectionAngle), ...
        0.5 * diameter * sin(sectionAngle), [0.8 0.8 0.8], ...
'EdgeColor','k','LineWidth', 1.2);
    title('Y = 0 longitudinal-plane velocity magnitude and streamlines');
    xlabel('X (m)');
    ylabel('Z (m)');
    xlim([-diameter, diameter]);
    ylim([-diameter, diameter]);

    %% Stage 3: Plot the z = 0 horizontal plane

    subplot(1, 2, 2);
    hold on;
    grid on;
    axis equal;

    middleZIndex = round(size(zGrid, 3) / 2);
    horizontalX = squeeze(xGrid(:, :, middleZIndex)); % [m]
    horizontalY = squeeze(yGrid(:, :, middleZIndex)); % [m]
    horizontalSpeed = squeeze(velocityField.mag(:, :, middleZIndex)); % [m/s]
    horizontalVelocityX = squeeze(velocityField.u(:, :, middleZIndex)); % [m/s]
    horizontalVelocityY = squeeze(velocityField.v(:, :, middleZIndex)); % [m/s]

    pcolor(horizontalX, horizontalY, horizontalSpeed);
    shading interp;
    colormap(gca, jet);
    colorbarHandle = colorbar;
    colorbarHandle.Label.String ='|V| / U_{\infty}';

    [streamlineStartX, streamlineStartY] = meshgrid(...
        min(horizontalX(:)), ...
        linspace(min(horizontalY(:)), max(horizontalY(:)), 15));
    streamline(horizontalX', horizontalY', horizontalVelocityX', ...
        horizontalVelocityY', streamlineStartX, streamlineStartY);
    fill(0.5 * diameter * cos(sectionAngle), ...
        0.5 * diameter * sin(sectionAngle), [0.8 0.8 0.8], ...
'EdgeColor','k','LineWidth', 1.2);

    title('Z = 0 horizontal-plane velocity magnitude and streamlines');
    xlabel('X (m)');
    ylabel('Y (m)');
    xlim([-diameter, diameter]);
    ylim([-diameter, diameter]);

    fprintf('[OK] Flow-field slice plots completed.\n');
end
