function mesh = read_bmf(filename)
% READ_BMF Execute the documented read_bmf operation.
%
% Syntax:
%   mesh = read_bmf(filename)
%
% Inputs:
%   filename        : [char|string] Input or output file path.
%
% Outputs:
%   mesh            : [struct] Boundary mesh with geometry expressed in SI units.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
%
%

    fid = fopen(filename, 'r');
    if fid == -1
        error('Unable to open BMF file: %s', filename);
    end
    file_cleanup = onCleanup(@() fclose(fid));

    header_line = fgetl(fid);

    line2 = fgetl(fid);
    c = textscan(line2, '%f %d');
    ulen = double(c{1});
    panel_type_flag = double(c{2});

    line3 = fgetl(fid);
    c = textscan(line3, '%d %d');
    isx = double(c{1});
    isy = double(c{2});

    line4 = fgetl(fid);
    c = textscan(line4, '%d');
    n_panels = double(c{1});

    raw_data = textscan(fid, '%f %f %f', 4 * n_panels);
    clear file_cleanup;  % closes fid via onCleanup

    if any(cellfun(@numel, raw_data) ~= 4 * n_panels)
        error('BMF vertex data are incomplete: %s', filename);
    end

    X = reshape(raw_data{1}, [4, n_panels])';
    Y = reshape(raw_data{2}, [4, n_panels])';
    Z = reshape(raw_data{3}, [4, n_panels])';

    vertices = zeros(n_panels, 4, 3);
    for j = 1:4
        vertices(:, j, 1) = X(:, j);
        vertices(:, j, 2) = Y(:, j);
        vertices(:, j, 3) = Z(:, j);
    end

    centers = zeros(n_panels, 3);
    normals = zeros(n_panels, 3);
    areas   = zeros(n_panels, 1);
    e1      = zeros(n_panels, 3);
    e2      = zeros(n_panels, 3);

    for k = 1:n_panels
        P = squeeze(vertices(k, :, :));
        centers(k, :) = mean(P, 1);

        % Area vector based on the two quadrilateral diagonals.
        d13 = P(3, :) - P(1, :);
        d24 = P(4, :) - P(2, :);
        n_raw = cross(d13, d24);
        mag = norm(n_raw);

        if mag < 1e-12
            n_raw = cross(P(2, :) - P(1, :), P(4, :) - P(1, :));
            mag = norm(n_raw);
        end
        if mag < 1e-12
            error('Degenerate panel %d in BMF file: %s', k, filename);
        end

        n = n_raw / mag;
        normals(k, :) = n;
        areas(k) = 0.5 * mag;

        % Project an edge onto the panel plane so that e1,e2,n form a
        % genuinely orthonormal right-handed coordinate system.
        candidates = [P(2, :) - P(1, :); ...
                      P(4, :) - P(1, :); ...
                      P(3, :) - P(1, :)];
        tangent = [0.0, 0.0, 0.0];
        for q = 1:size(candidates, 1)
            trial = candidates(q, :) - dot(candidates(q, :), n) * n;
            if norm(trial) > 1e-12
                tangent = trial;
                break;
            end
        end
        if norm(tangent) <= 1e-12
            error('Cannot construct local basis for panel %d: %s', k, filename);
        end

        e1(k, :) = tangent / norm(tangent);
        e2_vec = cross(n, e1(k, :));
        e2(k, :) = e2_vec / norm(e2_vec);
    end

    hydro = struct('Vx', 0.0, 'Vy', 0.0, 'Vz', 0.0, ...
                   'V_mean', 0.0, 'center_of_buoyancy', [0.0, 0.0, 0.0]);

    if panel_type_flag == 1
        sym_factor = double((1 + isx) * (1 + isy));

        hydro.Vx = sum(centers(:, 1) .* normals(:, 1) .* areas) * sym_factor;
        hydro.Vy = sum(centers(:, 2) .* normals(:, 2) .* areas) * sym_factor;
        hydro.Vz = sum(centers(:, 3) .* normals(:, 3) .* areas) * sym_factor;
        hydro.V_mean = (hydro.Vx + hydro.Vy + hydro.Vz) / 3.0;

        r_dot_n = sum(centers .* normals, 2);
        dV = (1.0 / 3.0) * r_dot_n .* areas;
        V_tot = sum(dV);

        if abs(V_tot) > 1e-12
            xb = (3.0 / 4.0) * sum(dV .* centers(:, 1)) / V_tot;
            yb = (3.0 / 4.0) * sum(dV .* centers(:, 2)) / V_tot;
            zb = (3.0 / 4.0) * sum(dV .* centers(:, 3)) / V_tot;
            if isx == 1, xb = 0.0; end
            if isy == 1, yb = 0.0; end
            hydro.center_of_buoyancy = [xb, yb, zb];
        end
    end

    mesh.header       = strtrim(header_line);
    mesh.ulen         = ulen;
    mesh.panel_type   = repmat(panel_type_flag, [n_panels, 1]);
    mesh.isx          = isx;
    mesh.isy          = isy;
    mesh.n_panels     = n_panels;
    mesh.vertices     = vertices;
    mesh.centers      = centers;
    mesh.normals      = normals;
    mesh.areas        = areas;
    mesh.e1           = e1;
    mesh.e2           = e2;
    mesh.hydrostatics = hydro;
end
