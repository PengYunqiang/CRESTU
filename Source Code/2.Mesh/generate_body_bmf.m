function mesh_body = generate_body_bmf(filename, diameter, isx, isy, panelDivisionCount)
% GENERATE_BODY_BMF Generate body bmf for the CRESTU hydrodynamic workflow.
%
% Syntax:
%   mesh_body = generate_body_bmf(filename, diameter, isx, isy, panelDivisionCount)
%
% Description:
%   Prepares boundary-panel geometry for the Rankine solver.
%   Global coordinates and panel-normal signs are preserved.
%
% Inputs:
%   filename           - [character vector or string scalar] Input or output file path.
%   diameter           - [scalar] Characteristic body diameter, [m].
%   isx                - [logical scalar] Symmetry-reduction flag for the x = 0 plane.
%   isy                - [logical scalar] Symmetry-reduction flag for the y = 0 plane.
%   panelDivisionCount - [scalar] Panel divisions per cube-face direction [-].
%
% Outputs:
%   mesh_body          - [struct] Generated body wetted-surface mesh and hydrostatics in SI units.
%
% Governing Equations / Theory:
%   Planar panel geometry, polygon moments, reflection transformations, and mesh-topology relations as applicable.
%
% References:
%   - Hess and Smith (1964); CRESTU BMF mesh-format and geometry conventions.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    if nargin < 1
        filename ='hemi_body.bmf';
    end
    if nargin < 2
        diameter = 10.0;
    end
    if nargin < 3
        isx = 1;
    end
    if nargin < 4
        isy = 1;
    end
    if nargin < 5
        panelDivisionCount = 8;
    end

    sphereRadius = diameter / 2.0; % [m]

    faces = {};
    [X1, Y1] = meshgrid(linspace(-1, 1, 2 * panelDivisionCount + 1), ...
        linspace(-1, 1, 2 * panelDivisionCount + 1));
    faces{end + 1} = {X1, Y1, -ones(size(X1))};
    [Y2, Z2] = meshgrid(linspace(-1, 1, 2 * panelDivisionCount + 1), ...
        linspace(-1, 0, panelDivisionCount + 1));
    faces{end + 1} = {ones(size(Y2)), Y2, Z2}; % +X
    faces{end + 1} = {-ones(size(Y2)), Y2, Z2}; % -X
    [X4, Z4] = meshgrid(linspace(-1, 1, 2 * panelDivisionCount + 1), ...
        linspace(-1, 0, panelDivisionCount + 1));
    faces{end + 1} = {X4, ones(size(X4)), Z4}; % +Y
    faces{end + 1} = {X4, -ones(size(X4)), Z4}; % -Y

    raw_panels = [];
    for faceIndex = 1:length(faces)
        Xf = faces{faceIndex}{1};
        Yf = faces{faceIndex}{2};
        Zf = faces{faceIndex}{3};
        [nr, nc] = size(Xf);
        for rowIndex = 1:nr - 1
            for columnIndex = 1:nc - 1
                c1 = [Xf(rowIndex, columnIndex), Yf(rowIndex, columnIndex), Zf(rowIndex, columnIndex)];
                c2 = [Xf(rowIndex, columnIndex + 1), Yf(rowIndex, columnIndex + 1), Zf(rowIndex, columnIndex + 1)];
                c3 = [Xf(rowIndex + 1, columnIndex + 1), Yf(rowIndex + 1, columnIndex + 1), Zf(rowIndex + 1, columnIndex + 1)];
                c4 = [Xf(rowIndex + 1, columnIndex), Yf(rowIndex + 1, columnIndex), Zf(rowIndex + 1, columnIndex)];

                p1 = project_sphere(c1, sphereRadius);
                p2 = project_sphere(c2, sphereRadius);
                p3 = project_sphere(c3, sphereRadius);
                p4 = project_sphere(c4, sphereRadius);

                p_center = (p1 + p2 + p3 + p4) / 4.0;
                if isx == 1 && (p_center(1) < -1e-6)
                    continue;
                end
                if isy == 1 && (p_center(2) < -1e-6)
                    continue;
                end

                v12 = p2 - p1;
                v14 = p4 - p1;
                if dot(cross(v12, v14), p_center) < 0
                    p_temp = p2;
                    p2 = p4;
                    p4 = p_temp;
                end
                raw_panels = cat(1, raw_panels, reshape([p1; p2; p3; p4], [1, 4, 3]));
            end
        end
    end

    mesh_body.header = sprintf('Pure-Quad Hemisphere D=%.2fm, ISX=%d, ISY=%d', diameter, isx, isy);
    mesh_body.ulen = 1.0;
    mesh_body.panel_type = ones(size(raw_panels, 1), 1); % Type 1 denotes the body boundary.
    mesh_body.isx = isx;
    mesh_body.isy = isy;
    mesh_body.n_panels = size(raw_panels, 1);
    mesh_body.vertices = raw_panels;

    write_bmf(filename, mesh_body);
    mesh_body = read_bmf(filename);
end

function projectedPoint = project_sphere(cubePoint, sphereRadius)
% PROJECT_SPHERE Project a Cartesian point radially onto a spherical surface.
%
% Syntax:
%   projectedPoint = project_sphere(cubePoint, sphereRadius)
%
% Description:
%   Prepares boundary-panel geometry for the Rankine solver.
%   Global coordinates and panel-normal signs are preserved.
%
% Inputs:
%   cubePoint          - [1 x 3] Cartesian cube-face point [-].
%   sphereRadius       - [scalar] Target sphere radius [m].
%
% Outputs:
%   projectedPoint     - [1 x 3] Projected Cartesian point [m].
%
% Governing Equations / Theory:
%   Planar panel geometry, polygon moments, reflection transformations, and mesh-topology relations as applicable.
%
% References:
%   - Hess and Smith (1964); CRESTU BMF mesh-format and geometry conventions.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    mappedX = tan(cubePoint(1) * pi / 4);
    mappedY = tan(cubePoint(2) * pi / 4);
    mappedZ = tan(cubePoint(3) * pi / 4);
    mappedRadius = sqrt(mappedX^2 + mappedY^2 + mappedZ^2);
    projectedPoint = sphereRadius * [mappedX, mappedY, mappedZ] / mappedRadius; % [m]
end
