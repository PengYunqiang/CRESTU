function mesh_body = generate_body_bmf(filename, diameter, isx, isy, N)
% GENERATE_BODY_BMF Execute the documented generate_body_bmf operation.
%
% Syntax:
%   mesh_body = generate_body_bmf(filename, diameter, isx, isy, N)
%
% Inputs:
%   filename        : [char|string] Input or output file path.
%   diameter        : [documented value] Input required by the implemented function contract.
%   isx             : [logical scalar] Reflection-symmetry flag for the x = 0 plane.
%   isy             : [logical scalar] Reflection-symmetry flag for the y = 0 plane.
%   N               : [integer scalar or array] Discrete count or index required by the algorithm.
%
% Outputs:
%   mesh_body       : [struct] Body boundary mesh with geometry expressed in SI units.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================

    if nargin < 1, filename = 'hemi_body.bmf'; end
    if nargin < 2, diameter = 10.0; end
    if nargin < 3, isx = 1; end
    if nargin < 4, isy = 1; end
    if nargin < 5, N = 8; end

    R = diameter / 2.0;

    faces = {};
    [X1, Y1] = meshgrid(linspace(-1, 1, 2*N+1), linspace(-1, 1, 2*N+1));
    faces{end+1} = {X1, Y1, -ones(size(X1))};
    [Y2, Z2] = meshgrid(linspace(-1, 1, 2*N+1), linspace(-1, 0, N+1));
    faces{end+1} = {ones(size(Y2)), Y2, Z2};  % +X
    faces{end+1} = {-ones(size(Y2)), Y2, Z2}; % -X
    [X4, Z4] = meshgrid(linspace(-1, 1, 2*N+1), linspace(-1, 0, N+1));
    faces{end+1} = {X4, ones(size(X4)), Z4};  % +Y
    faces{end+1} = {X4, -ones(size(X4)), Z4}; % -Y

    raw_panels = [];
    for f = 1:length(faces)
        Xf = faces{f}{1}; Yf = faces{f}{2}; Zf = faces{f}{3};
        [nr, nc] = size(Xf);
        for r = 1:nr-1
            for c = 1:nc-1
                c1 = [Xf(r, c),     Yf(r, c),     Zf(r, c)];
                c2 = [Xf(r, c+1),   Yf(r, c+1),   Zf(r, c+1)];
                c3 = [Xf(r+1, c+1), Yf(r+1, c+1), Zf(r+1, c+1)];
                c4 = [Xf(r+1, c),   Yf(r+1, c),   Zf(r+1, c)];

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

    mesh_body.header     = sprintf('Pure-Quad Hemisphere D=%.2fm, ISX=%d, ISY=%d', diameter, isx, isy);
    mesh_body.ulen       = 1.0;
    mesh_body.panel_type = ones(size(raw_panels, 1), 1); % Type 1 denotes the body boundary.
    mesh_body.isx        = isx;
    mesh_body.isy        = isy;
    mesh_body.n_panels   = size(raw_panels, 1);
    mesh_body.vertices   = raw_panels;

    write_bmf(filename, mesh_body);
    mesh_body = read_bmf(filename);
end

function P = project_sphere(pt, R)
% PROJECT_SPHERE Execute the documented project_sphere operation.
%
% Syntax:
%   P = project_sphere(pt, R)
%
% Inputs:
%   pt              : [documented value] Input required by the implemented function contract.
%   R               : [documented value] Input required by the implemented function contract.
%
% Outputs:
%   P               : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
    x = tan(pt(1) * pi / 4);
    y = tan(pt(2) * pi / 4);
    z = tan(pt(3) * pi / 4);
    P = R * [x, y, z] / sqrt(x^2 + y^2 + z^2);
end
