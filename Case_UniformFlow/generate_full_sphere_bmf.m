function mesh_body = generate_full_sphere_bmf(filename, diameter, isx, isy, panelDivisionCount)
% GENERATE_FULL_SPHERE_BMF Generate a closed cubed-sphere panel mesh and write a BMF file.
%
% Syntax:
%   mesh_body = generate_full_sphere_bmf(filename, diameter, isx, isy, panelDivisionCount)
%
% Description:
%   Generates a closed quadrilateral sphere by an equiangular cubed-sphere
%   projection. Panel vertices are ordered so normals point into the fluid.
%
% Inputs:
%   filename           - Output BMF path, character vector or string scalar [-].
%   diameter           - Sphere diameter, positive scalar [m].
%   isx                - Reflection flag for the x = 0 plane, 0 or 1 [-].
%   isy                - Reflection flag for the y = 0 plane, 0 or 1 [-].
%   panelDivisionCount - Panel divisions per cube-face direction [-].
%
% Outputs:
%   mesh_body          - Closed body-panel mesh with coordinates [m].
%
% Governing Equations / Theory:
%   Each cube coordinate is mapped with tan(pi*x/4), normalized, and
%   multiplied by the sphere radius. This removes pole-degenerate panels.
%
% References:
%   - Standard cubed-sphere mesh-generation methods.
%   - CRESTU theory and technical manual.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

    %% Stage 1: Apply defaults and validate inputs

    if nargin < 1 || isempty(filename)
        filename ='full_sphere_body.bmf';
    end

    if nargin < 2
        diameter = 10.0;
    end

    if nargin < 3
        isx = 0;
    end

    if nargin < 4
        isy = 0;
    end

    if nargin < 5
        panelDivisionCount = 8;
    end

    validateattributes(filename, {'char','string'}, {'nonempty'}, mfilename,'filename');
    validateattributes(diameter, {'double'}, {'scalar','positive','finite'}, mfilename,'diameter');
    validateattributes(isx, {'double','logical'}, {'scalar'}, mfilename,'isx');
    validateattributes(isy, {'double','logical'}, {'scalar'}, mfilename,'isy');
    validateattributes(panelDivisionCount, {'double'}, ...
        {'scalar','integer','positive','finite'}, mfilename,'panelDivisionCount');
    assert(ismember(isx, [0, 1]) && ismember(isy, [0, 1]), ...
'CRESTU:SphereSymmetryFlag','ISX and ISY must be 0 or 1.');

    sphereRadius = diameter / 2.0; % [m]

    %% Stage 2: Build the six cube faces

    gridCoordinates = linspace(-1, 1, 2 * panelDivisionCount + 1);
    [cubeCoordinateOne, cubeCoordinateTwo] = meshgrid(gridCoordinates, gridCoordinates);
    unitFace = ones(size(cubeCoordinateOne));

    faces = cell(6, 1);
    faces{1} = { unitFace, cubeCoordinateOne, cubeCoordinateTwo}; % +X face
    faces{2} = {-unitFace, -cubeCoordinateOne, cubeCoordinateTwo}; % -X face
    faces{3} = { cubeCoordinateOne, unitFace, cubeCoordinateTwo}; % +Y face
    faces{4} = {-cubeCoordinateOne, -unitFace, cubeCoordinateTwo}; % -Y face
    faces{5} = { cubeCoordinateOne, cubeCoordinateTwo, unitFace}; % +Z face
    faces{6} = {-cubeCoordinateOne, cubeCoordinateTwo, -unitFace}; % -Z face

    maximumPanelCount = 24 * panelDivisionCount ^ 2;
    rawPanels = zeros(maximumPanelCount, 4, 3); % [m]
    acceptedPanelCount = 0;

    %% Stage 3: Project all cube faces onto the sphere

    for faceIndex = 1:6
        faceX = faces{faceIndex}{1};
        faceY = faces{faceIndex}{2};
        faceZ = faces{faceIndex}{3};
        [faceRowCount, faceColumnCount] = size(faceX);

        for rowIndex = 1:faceRowCount - 1
            for columnIndex = 1:faceColumnCount - 1
                cubeVertexOne = [faceX(rowIndex, columnIndex), ...
                    faceY(rowIndex, columnIndex), faceZ(rowIndex, columnIndex)];
                cubeVertexTwo = [faceX(rowIndex, columnIndex + 1), ...
                    faceY(rowIndex, columnIndex + 1), faceZ(rowIndex, columnIndex + 1)];
                cubeVertexThree = [faceX(rowIndex + 1, columnIndex + 1), ...
                    faceY(rowIndex + 1, columnIndex + 1), faceZ(rowIndex + 1, columnIndex + 1)];
                cubeVertexFour = [faceX(rowIndex + 1, columnIndex), ...
                    faceY(rowIndex + 1, columnIndex), faceZ(rowIndex + 1, columnIndex)];

                vertexOne = project_sphere(cubeVertexOne, sphereRadius); % [m]
                vertexTwo = project_sphere(cubeVertexTwo, sphereRadius); % [m]
                vertexThree = project_sphere(cubeVertexThree, sphereRadius); % [m]
                vertexFour = project_sphere(cubeVertexFour, sphereRadius); % [m]
                panelCenter = (vertexOne + vertexTwo + vertexThree + vertexFour) / 4.0; % [m]

                if isx == 1 && panelCenter(1) < -1e-6
                    continue;
                end

                if isy == 1 && panelCenter(2) < -1e-6
                    continue;
                end

                edgeOne = vertexTwo - vertexOne; % [m]
                edgeTwo = vertexFour - vertexOne; % [m]

                if dot(cross(edgeOne, edgeTwo), panelCenter) < 0
                    temporaryVertex = vertexTwo;
                    vertexTwo = vertexFour;
                    vertexFour = temporaryVertex;
                end

                acceptedPanelCount = acceptedPanelCount + 1;
                rawPanels(acceptedPanelCount, :, :) = ...
                    reshape([vertexOne; vertexTwo; vertexThree; vertexFour], [1, 4, 3]);
            end
        end
    end

    rawPanels = rawPanels(1:acceptedPanelCount, :, :);

    %% Stage 4: Write the BMF mesh and check displaced volume

    totalPanelCount = size(rawPanels, 1);
    mesh_body = struct();
    mesh_body.header = sprintf('Pure-Quad Full Sphere D=%.2fm, ISX=%d, ISY=%d', ...
        diameter, isx, isy);
    mesh_body.ulen = 1.0;
    mesh_body.panel_type = ones(totalPanelCount, 1); % 1 = body panel
    mesh_body.isx = isx;
    mesh_body.isy = isy;
    mesh_body.n_panels = totalPanelCount;
    mesh_body.vertices = rawPanels;

    write_bmf(filename, mesh_body);
    mesh_body = read_bmf(filename);

    referenceVolume = 4.0 / 3.0 * pi * sphereRadius ^ 3; % [m^3]

    if isx == 1
        referenceVolume = referenceVolume / 2.0;
    end

    if isy == 1
        referenceVolume = referenceVolume / 2.0;
    end

    volumeErrorPercent = abs(mesh_body.hydrostatics.V_mean - referenceVolume) / ...
        referenceVolume * 100; % [%]
    fprintf('[OK] Full-sphere mesh generated | file = %s | panels = %d\n', ...
        filename, totalPanelCount);
    fprintf(['[OK] Volume check | reference = %.4f m^3 | integrated = %.4f m^3 | ', ...
'error = %.4f%%\n'], referenceVolume, mesh_body.hydrostatics.V_mean, ...
        volumeErrorPercent);
end

function projectedPoint = project_sphere(cubePoint, sphereRadius)
% PROJECT_SPHERE Project a cube-face point onto a sphere.
%
% Syntax:
%   projectedPoint = project_sphere(cubePoint, sphereRadius)
%
% Description:
%   Applies the equiangular tangent mapping and radial normalization.
%
% Inputs:
%   cubePoint    - Cube-face point, [1 x 3] [-].
%   sphereRadius - Sphere radius [m].
%
% Outputs:
%   projectedPoint - Projected sphere point, [1 x 3] [m].
%
% Governing Equations / Theory:
%   P = R*tan(pi*p/4)/norm(tan(pi*p/4)).
%
% References:
%   - Standard cubed-sphere mesh-generation methods.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

    %% Stage 1: Validate the projection inputs

    validateattributes(cubePoint, {'double'}, {'row','numel', 3,'finite'});
    validateattributes(sphereRadius, {'double'}, {'scalar','positive','finite'});

    %% Stage 2: Apply the equiangular tangent projection

    mappedX = tan(cubePoint(1) * pi / 4);
    mappedY = tan(cubePoint(2) * pi / 4);
    mappedZ = tan(cubePoint(3) * pi / 4);
    mappedRadius = sqrt(mappedX ^ 2 + mappedY ^ 2 + mappedZ ^ 2);
    projectedPoint = sphereRadius * [mappedX, mappedY, mappedZ] / mappedRadius; % [m]
end
