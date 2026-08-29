% UNIFORMFLOWMAIN Run the closed-sphere Hess-Smith uniform-flow benchmark.
%
% Description:
%   Generates a full-sphere panel mesh, solves the source strengths, checks
%   pressure and force conservation, exports a reconciliation CSV, and plots
%   surface and external flow fields. All physical inputs use SI units.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

clear;
clc;
close all;

%% Stage 1: Validate project dependencies and output paths

scriptDirectory = fileparts(mfilename('fullpath'));
projectDirectory = fileparts(scriptDirectory);
sourceDirectory = fullfile(projectDirectory,'Source Code');
assert(isfolder(sourceDirectory),'CRESTU:MissingSourceDirectory', ...
'Source directory was not found: %s', sourceDirectory);
addpath(genpath(sourceDirectory));
addpath(scriptDirectory);

requiredFunctions = {'generate_full_sphere_bmf','solve_uniform_flow', ...
'compute_flow_field','plot_field_slices','hess_smith_panel_velocity'};

for dependencyIndex = 1:numel(requiredFunctions)
    assert(exist(requiredFunctions{dependencyIndex},'file') == 2, ...
'CRESTU:MissingUniformFlowDependency', ...
'Required function was not found: %s', requiredFunctions{dependencyIndex});
end

meshFile = fullfile(scriptDirectory,'full_sphere_test.bmf');
summaryFile = fullfile(scriptDirectory,'UniformFlow_Reconciliation.csv');

fprintf('[INFO] Uniform-flow benchmark started.\n');

%% Stage 2: Generate the closed full-sphere body mesh

sphereDiameter = 10.0; % [m]
panelDivisionCount = 8; % [-], gives 6 * 4 * N^2 = 1536 panels
bodyMesh = generate_full_sphere_bmf(meshFile, sphereDiameter, 0, 0, ...
    panelDivisionCount);
assert(bodyMesh.n_panels > 0,'CRESTU:EmptyUniformFlowMesh', ...
'Generated body mesh contains no panels.');

%% Stage 3: Solve the Hess-Smith uniform-flow system

freeStreamVelocity = [1.0, 0.0, 0.0]; % [m/s]
flowResult = solve_uniform_flow(bodyMesh, freeStreamVelocity);
assert(numel(flowResult.Cp) == bodyMesh.n_panels, ...
'CRESTU:UniformFlowResultShape', ...
'Pressure coefficient must contain one value per panel.');

%% Stage 4: Check force balance and export reconciliation data

fluidDensity = 1025.0; % [kg/m^3]
dynamicPressure = 0.5 * fluidDensity * dot(freeStreamVelocity, ...
    freeStreamVelocity) .* flowResult.Cp; % [Pa]
hydrodynamicForce = zeros(1, 3); % [N]

for i = 1:bodyMesh.n_panels
    hydrodynamicForce = hydrodynamicForce - dynamicPressure(i) * ...
        bodyMesh.normals(i, :) * bodyMesh.areas(i);
end

maximumPressureCoefficient = max(flowResult.Cp); % [-]
minimumPressureCoefficient = min(flowResult.Cp); % [-]
reconciliation = table(bodyMesh.n_panels, sphereDiameter, ...
    norm(freeStreamVelocity), maximumPressureCoefficient, ...
    minimumPressureCoefficient, hydrodynamicForce(1), ...
    hydrodynamicForce(2), hydrodynamicForce(3), ...
'VariableNames', {'PanelCount','Diameter_m','FreeStreamSpeed_mps', ...
'CpMaximum','CpMinimum','ForceX_N','ForceY_N','ForceZ_N'});

fprintf('[INFO] Uniform-flow force and conservation summary:\n');
disp(reconciliation);
writetable(reconciliation, summaryFile);
fprintf('[OK] Reconciliation CSV exported | file = %s\n', summaryFile);

%% Stage 5: Plot the meridian pressure-coefficient comparison

panelCenters = bodyMesh.centers; % [m]
sliceTolerance = sphereDiameter / (2 * panelDivisionCount); % [m]
slicePanelIndices = find(abs(panelCenters(:, 2)) < sliceTolerance);
sliceX = panelCenters(slicePanelIndices, 1); % [m]
sliceZ = panelCenters(slicePanelIndices, 3); % [m]
slicePressureCoefficient = flowResult.Cp(slicePanelIndices); % [-]

polarAngleDegrees = atan2d(-sliceZ, -sliceX); % [deg]
[sortedPolarAngleDegrees, sortIndices] = sort(polarAngleDegrees);
sortedPressureCoefficient = slicePressureCoefficient(sortIndices);

figure('Color','w','Position', [100, 150, 780, 480],'Name','Cp Validation');
plot(sortedPolarAngleDegrees, sortedPressureCoefficient,'ro-', ...
'LineWidth', 1.5,'MarkerSize', 5.5, ...
'DisplayName','BEM solution (Hess-Smith)');
hold on;
grid on;

referencePolarAngleDegrees = linspace(-180, 180, 360); % [deg]
referencePressureCoefficient = 1.0 - ...
    2.25 * sind(referencePolarAngleDegrees) .^ 2; % [-]
plot(referencePolarAngleDegrees, referencePressureCoefficient,'b--', ...
'LineWidth', 2.0,'DisplayName', ...
'Analytical solution: C_p = 1 - 2.25 sin^2\theta');

xlabel(['Meridian polar angle \theta (deg) [0 deg: upstream, ', ...
'\pm90 deg: equator, \pm180 deg: downstream]'], ...
'FontSize', 11,'FontWeight','bold');
ylabel('Surface pressure coefficient C_p','FontSize', 11,'FontWeight','bold');
title('Closed-sphere surface C_p validation in uniform flow','FontSize', 12);
legend('Location','south');
ylim([-1.4, 1.2]);

%% Stage 6: Plot body pressure and surface velocity vectors

figureHandle = figure('Color','w','Position', [900, 150, 850, 650], ...
'Name','Flow 3D');
axesHandle = axes(figureHandle);
hold(axesHandle,'on');
grid(axesHandle,'on');
axis(axesHandle,'equal');
view(axesHandle, 135, 25);
xlabel(axesHandle,'X (m)');
ylabel(axesHandle,'Y (m)');
zlabel(axesHandle,'Z (m)');
title(axesHandle,'Body-surface C_p contours and velocity vectors');

surfaceX = squeeze(bodyMesh.vertices(:, :, 1))'; % [m]
surfaceY = squeeze(bodyMesh.vertices(:, :, 2))'; % [m]
surfaceZ = squeeze(bodyMesh.vertices(:, :, 3))'; % [m]
colormap(axesHandle, jet);
patch(axesHandle, surfaceX, surfaceY, surfaceZ, flowResult.Cp', ...
'FaceColor','flat','EdgeColor', [0.25 0.25 0.25],'LineWidth', 0.2);
colorbarHandle = colorbar(axesHandle);
colorbarHandle.Label.String ='C_p';

vectorStride = 1;
quiver3(axesHandle, bodyMesh.centers(1:vectorStride:end, 1), ...
    bodyMesh.centers(1:vectorStride:end, 2), ...
    bodyMesh.centers(1:vectorStride:end, 3), ...
    flowResult.V_total(1:vectorStride:end, 1), ...
    flowResult.V_total(1:vectorStride:end, 2), ...
    flowResult.V_total(1:vectorStride:end, 3), ...
    1.2,'k','LineWidth', 1.0);

%% Stage 7: Evaluate and plot external flow-field slices

domainHalfWidth = sphereDiameter; % [m]
gridPointCount = 41; % [-]
fprintf('[INFO] Generate the external flow grid | points per axis = %d\n', ...
    gridPointCount);
[xGrid, yGrid, zGrid, velocityField, pressureCoefficientField] = ...
    compute_flow_field(bodyMesh, flowResult.sigma, freeStreamVelocity, ...
    domainHalfWidth, gridPointCount);
plot_field_slices(xGrid, yGrid, zGrid, velocityField, ...
    pressureCoefficientField, bodyMesh, sphereDiameter);

fprintf('[OK] Uniform-flow benchmark completed.\n');
