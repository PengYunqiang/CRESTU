function mesh_box = generate_box_bmf(filename, L, B, D, isx, isy, Nx, Ny, Nz)
% GENERATE_BOX_BMF Execute the documented generate_box_bmf operation.
%
% Syntax:
%   mesh_box = generate_box_bmf(filename, L, B, D, isx, isy, Nx, Ny, Nz)
%
% Inputs:
%   filename        : [char|string] Input or output file path.
%   L               : [documented value] Input required by the implemented function contract.
%   B               : [documented value] Input required by the implemented function contract.
%   D               : [documented value] Input required by the implemented function contract.
%   isx             : [logical scalar] Reflection-symmetry flag for the x = 0 plane.
%   isy             : [logical scalar] Reflection-symmetry flag for the y = 0 plane.
%   Nx              : [integer scalar or array] Discrete count or index required by the algorithm.
%   Ny              : [integer scalar or array] Discrete count or index required by the algorithm.
%   Nz              : [integer scalar or array] Discrete count or index required by the algorithm.
%
% Outputs:
%   mesh_box        : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
%

    if nargin < 1 || isempty(filename), filename = 'box_body.bmf'; end
    if nargin < 2, L = 20.0; end
    if nargin < 3, B = 10.0; end
    if nargin < 4, D = 4.0;  end
    if nargin < 5, isx = 0;  end
    if nargin < 6, isy = 0;  end
    if nargin < 7, Nx = 10;  end
    if nargin < 8, Ny = 6;   end
    if nargin < 9, Nz = 5;   end

    if isx == 1
        x_range = [0, L/2];
    else
        x_range = [-L/2, L/2];
    end

    if isy == 1
        y_range = [0, B/2];
    else
        y_range = [-B/2, B/2];
    end

    z_range = [-D, 0];

    x_lin = linspace(x_range(1), x_range(2), Nx + 1);
    y_lin = linspace(y_range(1), y_range(2), Ny + 1);
    z_lin = linspace(z_range(1), z_range(2), Nz + 1);

    raw_panels = [];

    for i = 1:Nx
        for j = 1:Ny
            x1 = x_lin(i);   x2 = x_lin(i+1);
            y1 = y_lin(j);   y2 = y_lin(j+1);

            p1 = [x1, y1, -D];
            p2 = [x1, y2, -D];
            p3 = [x2, y2, -D];
            p4 = [x2, y1, -D];

            raw_panels = cat(1, raw_panels, reshape([p1; p2; p3; p4], [1, 4, 3]));
        end
    end

    for j = 1:Ny
        for k = 1:Nz
            y1 = y_lin(j);   y2 = y_lin(j+1);
            z1 = z_lin(k);   z2 = z_lin(k+1);

            p1 = [L/2, y1, z1];
            p2 = [L/2, y2, z1];
            p3 = [L/2, y2, z2];
            p4 = [L/2, y1, z2];

            raw_panels = cat(1, raw_panels, reshape([p1; p2; p3; p4], [1, 4, 3]));
        end
    end

    if isx == 0
        for j = 1:Ny
            for k = 1:Nz
                y1 = y_lin(j);   y2 = y_lin(j+1);
                z1 = z_lin(k);   z2 = z_lin(k+1);

                p1 = [-L/2, y2, z1];
                p2 = [-L/2, y1, z1];
                p3 = [-L/2, y1, z2];
                p4 = [-L/2, y2, z2];

                raw_panels = cat(1, raw_panels, reshape([p1; p2; p3; p4], [1, 4, 3]));
            end
        end
    end

    for i = 1:Nx
        for k = 1:Nz
            x1 = x_lin(i);   x2 = x_lin(i+1);
            z1 = z_lin(k);   z2 = z_lin(k+1);

            p1 = [x2, B/2, z1];
            p2 = [x1, B/2, z1];
            p3 = [x1, B/2, z2];
            p4 = [x2, B/2, z2];

            raw_panels = cat(1, raw_panels, reshape([p1; p2; p3; p4], [1, 4, 3]));
        end
    end

    if isy == 0
        for i = 1:Nx
            for k = 1:Nz
                x1 = x_lin(i);   x2 = x_lin(i+1);
                z1 = z_lin(k);   z2 = z_lin(k+1);

                p1 = [x1, -B/2, z1];
                p2 = [x2, -B/2, z1];
                p3 = [x2, -B/2, z2];
                p4 = [x1, -B/2, z2];

                raw_panels = cat(1, raw_panels, reshape([p1; p2; p3; p4], [1, 4, 3]));
            end
        end
    end

    total_panels = size(raw_panels, 1);
    mesh_box.header     = sprintf('Box Barge L=%.1fm, B=%.1fm, D=%.1fm (ISX=%d, ISY=%d)', L, B, D, isx, isy);
    mesh_box.ulen       = 1.0;
    mesh_box.panel_type = ones(total_panels, 1); % Type 1 denotes the body boundary.
    mesh_box.isx        = isx;
    mesh_box.isy        = isy;
    mesh_box.n_panels   = total_panels;
    mesh_box.vertices   = raw_panels;

    write_bmf(filename, mesh_box);

    mesh_box = read_bmf(filename);

    V_exact = L * B * D;
    Zb_exact = -D / 2.0;
    fprintf('\n========== 长方体浮箱静水力理论对标 ==========\n');
    fprintf(' 浮箱尺寸     : L = %.2f m, B = %.2f m, D = %.2f m\n', L, B, D);
    fprintf(' 理论排水体积 : %12.4f m^3\n', V_exact);
    fprintf(' 数值积分体积 : %12.4f m^3 (误差: %.4f%%)\n', mesh_box.hydrostatics.Vz, ...
        abs(mesh_box.hydrostatics.Vz - V_exact)/V_exact*100);
    fprintf(' 理论浮心 Zb  : %12.4f m\n', Zb_exact);
    fprintf(' 数值浮心 Zb  : %12.4f m (误差: %.4f%%)\n', ...
        mesh_box.hydrostatics.center_of_buoyancy(3), ...
        abs(mesh_box.hydrostatics.center_of_buoyancy(3) - Zb_exact)/abs(Zb_exact)*100);
    fprintf('==============================================\n\n');
end
