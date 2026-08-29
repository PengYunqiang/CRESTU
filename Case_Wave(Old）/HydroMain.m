clear;
clc;
close all;

% =========================================================================
% 1. Read main control parameters
% =========================================================================
cfg = read_config('CRESTU.cfg');

% =========================================================================
% 2. Assemble body boundaries and transform them to the global frame
% =========================================================================
body_list = cell(cfg.n_bodies, 1);
total_body_panels = 0;
for bodyIndex = 1:cfg.n_bodies
    b_cfg = cfg.bodies(bodyIndex);
    generate_body_bmf(b_cfg.mesh_file, 10.0, cfg.isx, cfg.isy, 5);
    mesh_loc = read_bmf(b_cfg.mesh_file);
    mesh_glob = transform_body_mesh(mesh_loc, b_cfg);
    body_list{bodyIndex} = mesh_glob;
    total_body_panels = total_body_panels + mesh_glob.n_panels;
end

% 3. Extract the global waterline
waterline = extract_waterline(body_list{1}, cfg.z_tol);

% 4. Generate the free-surface mesh and export BMF type 2
mesh_fs = generate_free_surface_bmf(cfg.files.fs, waterline, cfg);

% 5. Build finite-depth seabed type 4 and far-field type 5 boundaries
has_seabed = (cfg.water_depth > 0);
mesh_seabed = [];
mesh_farfield = [];
if has_seabed
% Generate the closed seabed mesh
    mesh_seabed = generate_seabed_mesh(waterline, cfg, mesh_fs);
    write_bmf(cfg.files.seabed, mesh_seabed);
    mesh_farfield = generate_farfield_mesh(waterline, cfg, cfg.fs.nz_farfield);
    write_bmf(cfg.files.farfield, mesh_farfield);
    total_dofs = total_body_panels + mesh_fs.n_panels + ...
        mesh_seabed.n_panels + mesh_farfield.n_panels;
else
    total_dofs = total_body_panels + mesh_fs.n_panels;
end

% 6. Assemble domain statistics
stats.total_body_panels = total_body_panels;
stats.fs_panels = mesh_fs.n_panels;
stats.seabed_panels = ternary(has_seabed, mesh_seabed.n_panels, 0);
stats.farfield_panels = ternary(has_seabed, mesh_farfield.n_panels, 0);
stats.total_dofs = total_dofs;

domain.cfg = cfg;
domain.body_list = body_list;
domain.fs = mesh_fs;
domain.seabed = mesh_seabed;
domain.farfield = mesh_farfield;
domain.waterline = waterline;
domain.stats = stats;

% Print the fluid-domain assembly report
fprintf('================= BMF domain assembly summary =================\n');
fprintf(' Case name        : %s\n', cfg.case_name);
fprintf(' Body count (NBODY): %d (total body panels: %d)\n', cfg.n_bodies, total_body_panels);
fprintf(' Free-surface panels  : %d (outer radius R_out = %.1f m)\n', mesh_fs.n_panels, cfg.fs.r_outer);
if has_seabed
    fprintf(' Seabed panels  : %d (water depth h = %.1f m)\n', mesh_seabed.n_panels, cfg.water_depth);
    fprintf(' Far-field panels  : %d\n', mesh_farfield.n_panels);
end
fprintf(' --------------------------------------------------------\n');
fprintf(' Reduced-domain unknowns: %d\n', total_dofs);
if ~isempty(body_list{1}.hydrostatics)
    fprintf(' Displaced volume Vz: %.4f m^3, center of buoyancy Zb: %.4f m\n', ...
        body_list{1}.hydrostatics.Vz, body_list{1}.hydrostatics.center_of_buoyancy(3));
end
fprintf('========================================================\n\n');

% Plot the complete three-dimensional fluid domain
plot_bmf_domain(domain,'wireframe');

% =========================================================================
% 7. Extract domain and body geometry arrays
% =========================================================================
centers = [];
normals = [];
verts = [];
body_areas = [];
for bodyIndex = 1:cfg.n_bodies
    centers = [centers; domain.body_list{bodyIndex}.centers];
    normals = [normals; domain.body_list{bodyIndex}.normals];
    verts = [verts; domain.body_list{bodyIndex}.vertices];
    body_areas = [body_areas; domain.body_list{bodyIndex}.areas];
end
centers = [centers; domain.fs.centers];
normals = [normals; domain.fs.normals];
verts = [verts; domain.fs.vertices];
if has_seabed
    centers = [centers; domain.seabed.centers];
    normals = [normals; domain.seabed.normals];
    verts = [verts; domain.seabed.vertices];

    centers = [centers; domain.farfield.centers];
    normals = [normals; domain.farfield.normals];
    verts = [verts; domain.farfield.vertices];
end

body_centers = centers(1:domain.stats.total_body_panels, :);
body_normals = normals(1:domain.stats.total_body_panels, :);
nj = compute_generalized_normals(body_centers, body_normals, domain.body_list);

% =========================================================================
% 8. Run the frequency-domain solver and load post-processing
% =========================================================================
n_freqs = length(cfg.freq.omegas);
wave_headings = cfg.wave.headings; % Wave headings from the validated configuration.
n_headings = length(wave_headings);

results(n_freqs) = struct();
pot_file_cache = sprintf('%s_Potential_Cache.mat', cfg.case_name);

for i_f = 1:n_freqs
    omega = cfg.freq.omegas(i_f);
    fprintf('\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n');
    fprintf(' [INFO] Start frequency case %d/%d: omega = %.4f rad/s\n', i_f, n_freqs, omega);
    fprintf('<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n');

    phi_rad_body = [];
    phi_diff_body = [];
    phi_I_body = [];

% --- 8.1 Potential solution controlled by IPOTEN ---
    if cfg.run.ipoten == 1
        [A, S_body] = assemble_rankine_matrix(total_dofs, centers, normals, verts, domain.stats, omega, cfg.grav, cfg);

        fprintf('[INFO] Factor the system matrix with LU decomposition (O(N^3))...\n');
        tic;
        [L, U, P] = lu(A);
        fprintf('[OK] LU factorization completed. Elapsed time: %.2f s\n', toc);

% Solve radiation potentials
        if cfg.run.irad == 1
            fprintf('[INFO] Solve all radiation modes...\n');
            b_rad = S_body * (1i * omega * nj);
            phi_rad = U \ (L \ (P * b_rad));
            phi_rad_body = phi_rad(1:domain.stats.total_body_panels, :);
        end

% Solve diffraction potentials
        if cfg.run.idiff == 1
            fprintf('[INFO] Solve all diffraction headings...\n');
            b_diff = complex(zeros(total_dofs, n_headings));
            phi_I_body = complex(zeros(domain.stats.total_body_panels, n_headings));
            for headingIndex = 1:n_headings
                [phi_I, dphi_I_dn] = compute_incident_wave(body_centers, body_normals, ...
                                        omega, cfg.grav, cfg.water_depth, wave_headings(headingIndex), 1.0);
                phi_I_body(:, headingIndex) = phi_I;
                b_diff(:, headingIndex) = S_body * (-dphi_I_dn);
            end
            phi_diff = U \ (L \ (P * b_diff));
            phi_diff_body = phi_diff(1:domain.stats.total_body_panels, :);
        end

% Save frequency potentials to the cache
% The cache may also be stored by frequency when required.
    else
        fprintf('[INFO] IPOTEN = 0. Load potentials from the local cache...\n');
        if exist(pot_file_cache,'file')
            load(pot_file_cache,'saved_pots');
            phi_rad_body = saved_pots(i_f).phi_rad_body;
            phi_diff_body = saved_pots(i_f).phi_diff_body;
            phi_I_body = saved_pots(i_f).phi_I_body;
        else
            error('Potential cache %s was not found. Run once with IPOTEN = 1.', pot_file_cache);
        end
    end

% --- 8.2 Load, added-mass, and damping post-processing controlled by IFORCE ---
    A_mat = [];
    B_mat = [];
    F_exc = [];
    if cfg.run.iforce == 1
        fprintf('[INFO] Compute hydrodynamic loads and coefficients...\n');

% Radiation added mass and damping
        if cfg.run.irad == 1 && ~isempty(phi_rad_body)
            [A_mat, B_mat] = compute_hydrodynamic_coeffs(phi_rad_body, nj, body_centers, body_areas, omega, cfg.grav, cfg.rho);
        end

% Diffraction excitation force
        if cfg.run.idiff == 1 && ~isempty(phi_diff_body)
            F_exc = complex(zeros(6 * cfg.n_bodies, n_headings));
            for headingIndex = 1:n_headings
                if isempty(phi_I_body)
                    [phi_I, ~] = compute_incident_wave(body_centers, body_normals, ...
                                    omega, cfg.grav, cfg.water_depth, wave_headings(headingIndex), 1.0);
                    phi_I_col = phi_I;
                else
                    phi_I_col = phi_I_body(:, headingIndex);
                end

                phi_total = phi_I_col + phi_diff_body(:, headingIndex);
                pressure = 1i * omega * cfg.rho * phi_total;
                for dof = 1:(6 * cfg.n_bodies)
                    F_exc(dof, headingIndex) = sum(pressure .* nj(:, dof) .* body_areas);
                end
            end
        end
    else
        fprintf('[WARN] IFORCE = 0. Skip load calculations.\n');
    end

% 8.3 Assemble results
    results(i_f).omega = omega;
    results(i_f).A_mat = A_mat;
    results(i_f).B_mat = B_mat;
    results(i_f).F_exc = F_exc;

    fprintf('[OK] Frequency %.4f rad/s completed.\n', omega);
end

fprintf('\n================== Frequency-domain analysis completed ==================\n');

% Local utility functions
function selectedValue = ternary(condition, trueValue, falseValue)
% TERNARY Select one of two values from a logical condition.
%
% Syntax:
%   out = ternary(cond, a, bodyIndex)
%
% Description:
%   Select one of two values from a logical condition.
%   The implementation preserves the CRESTU solver and data-field contracts.
%
% Inputs:
%   cond - Logical selection condition [-].
%   a - Value returned when the condition is true [-].
%   bodyIndex - Value returned when the condition is false [-].
%
% Outputs:
%   out - Selected input value [-].
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

    %% Stage 2: Select the requested value

    if condition
        selectedValue = trueValue;
    else
        selectedValue = falseValue;
    end
end
