function mesh = read_bmf(filename)
% READ_BMF Read bmf for the CRESTU hydrodynamic workflow.
%
% Syntax:
%   mesh = read_bmf(filename)
%
% Description:
%   Prepares boundary-panel geometry for the Rankine solver.
%   Global coordinates and panel-normal signs are preserved.
%
% Inputs:
%   filename           - [character vector or string scalar] Input or output file path.
%
% Outputs:
%   mesh               - [struct] Generated boundary-panel mesh with coordinates in [m] and areas in [m^2].
%
% Governing Equations / Theory:
%   Planar panel geometry, polygon moments, reflection transformations, and mesh-topology relations as applicable.
%
% References:
%   - Hess and Smith (1964); CRESTU BMF mesh-format and geometry conventions.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    fid = fopen(filename,'r');
    if fid == -1
        error('Unable to open BMF file: %s', filename);
    end
    file_cleanup = onCleanup(@() fclose(fid));

    header_line = fgetl(fid);

    line2 = fgetl(fid);
    parsedLine = textscan(line2,'%f %d');
    ulen = double(parsedLine{1});
    panel_type_flag = double(parsedLine{2});

    line3 = fgetl(fid);
    parsedLine = textscan(line3,'%d %d');
    isx = double(parsedLine{1});
    isy = double(parsedLine{2});

    line4 = fgetl(fid);
    parsedLine = textscan(line4,'%d');
    n_panels = double(parsedLine{1});

    raw_data = textscan(fid,'%f %f %f', 4 * n_panels);
    clear file_cleanup; % closes fid via onCleanup

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
    areas = zeros(n_panels, 1);
    e1 = zeros(n_panels, 3);
    e2 = zeros(n_panels, 3);

    for k = 1:n_panels
        panelVertices = squeeze(vertices(k, :, :));
        centers(k, :) = mean(panelVertices, 1);

% Area vector based on the two quadrilateral diagonals.
        d13 = panelVertices(3, :) - panelVertices(1, :);
        d24 = panelVertices(4, :) - panelVertices(2, :);
        n_raw = cross(d13, d24);
        mag = norm(n_raw);

        if mag < 1e-12
            n_raw = cross(panelVertices(2, :) - panelVertices(1, :), ...
                panelVertices(4, :) - panelVertices(1, :));
            mag = norm(n_raw);
        end
        if mag < 1e-12
            error('Degenerate panel %d in BMF file: %s', k, filename);
        end

        unitNormal = n_raw / mag;
        normals(k, :) = unitNormal;
        areas(k) = 0.5 * mag;

% Project an edge onto the panel plane so that e1,e2,unitNormal form a
% genuinely orthonormal right-handed coordinate system.
        candidates = [panelVertices(2, :) - panelVertices(1, :); ...
                      panelVertices(4, :) - panelVertices(1, :); ...
                      panelVertices(3, :) - panelVertices(1, :)];
        tangent = [0.0, 0.0, 0.0];
        for candidateIndex = 1:size(candidates, 1)
            trial = candidates(candidateIndex, :) - dot(candidates(candidateIndex, :), unitNormal) * unitNormal;
            if norm(trial) > 1e-12
                tangent = trial;
                break;
            end
        end
        if norm(tangent) <= 1e-12
            error('Cannot construct local basis for panel %d: %s', k, filename);
        end

        e1(k, :) = tangent / norm(tangent);
        e2_vec = cross(unitNormal, e1(k, :));
        e2(k, :) = e2_vec / norm(e2_vec);
    end

    hydro = struct('Vx', 0.0,'Vy', 0.0,'Vz', 0.0, ...
'V_mean', 0.0,'center_of_buoyancy', [0.0, 0.0, 0.0]);

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
            if isx == 1
                xb = 0.0;
            end
            if isy == 1
                yb = 0.0;
            end
            hydro.center_of_buoyancy = [xb, yb, zb];
        end
    end

    mesh.header = strtrim(header_line);
    mesh.ulen = ulen;
    mesh.panel_type = repmat(panel_type_flag, [n_panels, 1]);
    mesh.isx = isx;
    mesh.isy = isy;
    mesh.n_panels = n_panels;
    mesh.vertices = vertices;
    mesh.centers = centers;
    mesh.normals = normals;
    mesh.areas = areas;
    mesh.e1 = e1;
    mesh.e2 = e2;
    mesh.hydrostatics = hydro;
end
