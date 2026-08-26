function mesh_body = generate_body_bmf(filename, diameter, isx, isy, N)
% GENERATE_BODY_BMF Generate body bmf for the CRESTU hydrodynamic workflow.
%
% Syntax:
%   mesh_body = generate_body_bmf(filename, diameter, isx, isy, N)
%
% Description:
%   The routine constructs, transforms, validates, or visualizes boundary-panel geometry used by the Rankine solver. Coordinates are expressed in the global Cartesian frame and panel orientation is preserved so that normals remain consistent with boundary-integral signs.
%
% Inputs:
%   filename           - [character vector or string scalar] Input or output file path.
%   diameter           - [scalar] Characteristic body diameter, [m].
%   isx                - [logical scalar] Symmetry-reduction flag for the x = 0 plane.
%   isy                - [logical scalar] Symmetry-reduction flag for the y = 0 plane.
%   N                  - [scalar] Total number of boundary panels, dimensionless.
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

%% --- 1. Validate Inputs and Initialize the Algorithm ---

    if nargin < 1, filename = 'hemi_body.bmf'; end
    if nargin < 2, diameter = 10.0; end
    if nargin < 3, isx = 1; end
    if nargin < 4, isy = 1; end
    if nargin < 5, N = 8; end

    R = diameter / 2.0;

    faces = {};
    [X1, Y1] = meshgrid(linspace(-1, 1, 2 * N + 1), linspace(-1, 1, 2 * N + 1));
    faces{end + 1} = {X1, Y1, -ones(size(X1))};
    [Y2, Z2] = meshgrid(linspace(-1, 1, 2 * N + 1), linspace(-1, 0, N + 1));
    faces{end + 1} = {ones(size(Y2)), Y2, Z2}; % +X
    faces{end + 1} = {-ones(size(Y2)), Y2, Z2}; % -X
    [X4, Z4] = meshgrid(linspace(-1, 1, 2 * N + 1), linspace(-1, 0, N + 1));
    faces{end + 1} = {X4, ones(size(X4)), Z4}; % +Y
    faces{end + 1} = {X4, -ones(size(X4)), Z4}; % -Y

    raw_panels = [];
    for f = 1:length(faces)
        Xf = faces{f}{1}; Yf = faces{f}{2}; Zf = faces{f}{3};
        [nr, nc] = size(Xf);
        for r = 1:nr - 1
            for c = 1:nc - 1
                c1 = [Xf(r, c), Yf(r, c), Zf(r, c)];
                c2 = [Xf(r, c + 1), Yf(r, c + 1), Zf(r, c + 1)];
                c3 = [Xf(r + 1, c + 1), Yf(r + 1, c + 1), Zf(r + 1, c + 1)];
                c4 = [Xf(r + 1, c), Yf(r + 1, c), Zf(r + 1, c)];

                p1 = project_sphere(c1, R);
                p2 = project_sphere(c2, R);
                p3 = project_sphere(c3, R);
                p4 = project_sphere(c4, R);

                p_center = (p1 + p2 + p3 + p4) / 4.0;
                if isx == 1 && (p_center(1) < -1e-6), continue; end
                if isy == 1 && (p_center(2) < -1e-6), continue; end

                v12 = p2 - p1; v14 = p4 - p1;
                if dot(cross(v12, v14), p_center) < 0
                    p_temp = p2; p2 = p4; p4 = p_temp;
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

function P = project_sphere(pt, R)
% PROJECT_SPHERE Project a Cartesian point radially onto a spherical surface.
%
% Syntax:
%   P = project_sphere(pt, R)
%
% Description:
%   The routine constructs, transforms, validates, or visualizes boundary-panel geometry used by the Rankine solver. Coordinates are expressed in the global Cartesian frame and panel orientation is preserved so that normals remain consistent with boundary-integral signs.
%
% Inputs:
%   pt                 - [1 x 3] Cartesian point to project, [m].
%   R                  - [scalar] Target sphere radius, [m].
%
% Outputs:
%   P                  - [1 x 3] Projected Cartesian point, [m].
%
% Governing Equations / Theory:
%   Planar panel geometry, polygon moments, reflection transformations, and mesh-topology relations as applicable.
%
% References:
%   - Hess and Smith (1964); CRESTU BMF mesh-format and geometry conventions.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% --- 1. Validate Inputs and Initialize the Algorithm ---

    x = tan(pt(1) * pi / 4);
    y = tan(pt(2) * pi / 4);
    z = tan(pt(3) * pi / 4);
    P = R * [x, y, z] / sqrt(x^2 + y^2 + z^2);
end
