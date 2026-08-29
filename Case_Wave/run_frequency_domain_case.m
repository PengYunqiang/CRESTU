function results = run_frequency_domain_case(config_file)
% RUN_FREQUENCY_DOMAIN_CASE Solve radiation, diffraction, response, and mean-drift loads in the frequency domain.
%
% Syntax:
%   results = run_frequency_domain_case(config_file)
%
% Description:
%   Solve radiation, diffraction, response, and mean-drift loads in the frequency domain.
%   The implementation preserves the CRESTU solver and data-field contracts.
%
% Inputs:
%   config_file - Path to a CRESTU configuration file, character vector or string scalar [-].
%
% Outputs:
%   results - Complete CRESTU frequency-domain result structure [-].
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
    arguments
        config_file {mustBeTextScalar}
    end

    %% Stage 1: Validate the case input and assemble the base domain

    assert(isfile(config_file),'CRESTU:MissingConfiguration', ...
'Configuration file was not found: %s', config_file);
    fprintf('[INFO] Build the hydrodynamic domain | config = %s\n', config_file);
    domain = build_bmf_domain(config_file);
    cfg = domain.cfg;
    geom = domain.geometry;
    nb = domain.stats.total_body_panels;
    ndof = 6 * cfg.n_bodies;
    nf = numel(cfg.freq.omegas);
    nh = numel(cfg.wave.headings);
    assert(nb > 0,'CRESTU:EmptyBodyBoundary', ...
'At least one body panel is required.');
    assert(nf > 0 && all(cfg.freq.omegas > 0),'CRESTU:InvalidFrequencyGrid', ...
'All wave angular frequencies must be positive.');

    % A five-module model has 30 global DOFs. Module m uses
    % 6*(m-1)+(1:6). Local DOF 1-6: Surge, Sway, Heave, Roll, Pitch, Yaw.
    body_centers = geom.centers(1:nb, :); % [m]
    body_normals = geom.normals(1:nb, :); % [-]
    body_areas = geom.areas(1:nb); % [m^2]
    nj = compute_generalized_normals(body_centers, body_normals, domain.body_list);
    assert(isequal(size(nj), [nb, ndof]),'CRESTU:GeneralizedNormalShape', ...
'Generalized normals must have size [N_body_panels x 6*N_bodies].');

    %% Stage 2: Build symmetry groups and hydrostatic restoring terms

    mode_parity = get_mode_parities(cfg.n_bodies, cfg.isx, cfg.isy);
    group_parity = unique(mode_parity,'rows','stable');
    ng = size(group_parity, 1);
    group_modes = cell(ng, 1);

    for symmetryGroupIndex = 1:ng
        group_modes{symmetryGroupIndex} = find(all(...
            mode_parity == group_parity(symmetryGroupIndex, :), 2));
    end

    symmetry = struct('isx', cfg.isx,'isy', cfg.isy, ...
'multiplicity', cfg.symmetry.multiplicity, ...
'mode_parity', mode_parity,'group_parity', group_parity, ...
'reduced_dofs', geom.total_panels, ...
'full_equivalent_dofs', geom.total_panels * cfg.symmetry.multiplicity);
    [C, hydrostatics] = compute_hydrostatic_matrix(domain.body_list, cfg);
    assert(isequal(size(C), [ndof, ndof]),'CRESTU:HydrostaticMatrixShape', ...
'Hydrostatic matrix must have size [6*N_bodies x 6*N_bodies].');

    %% Stage 3: Initialize or load the potential cache

    if cfg.run.ipoten == 0
        pot_cache = load_potential_cache(cfg.files.potential_cache,cfg,geom);
        if pot_cache.schema_version ~= 4
            error('CRESTU:CacheSchema','Solver requires potential cache schema 4 (frequency-local outer meshes).');
        end
    else
        empty_entry = struct('omega',0,'phi_rad_body',complex(zeros(nb,ndof)), ...
'phi_diff_body',complex(zeros(nb,nh)),'phi_incident_body',complex(zeros(nb,nh)), ...
'phi_diff_components',complex(zeros(nb,nh,ng)), ...
'phi_incident_components',complex(zeros(nb,nh,ng)), ...
'frequency_domain_stats',struct(),'frequency_environment',zeros(1,3));
        pot_cache = struct('schema_version',4,'case_name',cfg.case_name,'omegas',cfg.freq.omegas, ...
'headings',cfg.wave.headings,'n_body_panels',nb,'n_total_panels',geom.total_panels, ...
'environment',environment_signature(cfg),'geometry_signature',geometry_signature(geom), ...
'mode_parity',mode_parity,'group_parity',group_parity,'n_dofs',ndof, ...
'created',char(datetime('now','Format','yyyyMMdd''T''HHmmss')), ...
'entries',repmat(empty_entry,nf,1));
    end

    A_all = zeros(ndof,ndof,nf);
    B_all = zeros(ndof,ndof,nf);
    F_all = complex(zeros(ndof,nh,nf));
    pressure_all = complex(zeros(nb,nh,nf));
    diagnostics = repmat(empty_diagnostics(ndof),nf,1);
    assembly_info = cell(nf,ng);

    %% Stage 4: Solve potentials and first-order loads by frequency

    for k = 1:nf
        omega = cfg.freq.omegas(k); % [rad/s]
        fprintf('[INFO] Frequency %d/%d | omega = %.6g rad/s\n', k, nf, omega);
        if cfg.run.ipoten == 1
            frequency_domain = build_frequency_local_domain(domain,omega);
            frequency_geometry = frequency_domain.geometry;
            frequency_cfg = frequency_domain.cfg;
            entry = pot_cache.entries(k);
            entry.omega = omega;
            entry.frequency_domain_stats = frequency_domain.stats;
            entry.frequency_environment = [frequency_cfg.fs.r_inner,frequency_cfg.fs.r_outer, ...
                frequency_cfg.fs.tuning_report.dynamic_mu0];
            for symmetryGroupIndex = 1:ng
                parity = group_parity(symmetryGroupIndex, :);
                modes = group_modes{symmetryGroupIndex};
                [K,S,assembly_info{k,symmetryGroupIndex}] = assemble_rankine_matrix(frequency_geometry.total_panels, ...
                    frequency_geometry.centers,frequency_geometry.normals,frequency_geometry.vertices, ...
                    frequency_domain.stats,omega,frequency_cfg.grav,frequency_cfg,parity);
                factors = [];
                if cfg.run.irad == 1
                    rhs = S * (1i * omega * nj(:,modes));
                    [phi,factors] = solve_complex_system(K,rhs,factors);
                    entry.phi_rad_body(:,modes) = phi(1:nb,:);
                end
                if cfg.run.idiff == 1
                    phi_I = complex(zeros(nb,nh));
                    rhs = complex(zeros(frequency_geometry.total_panels,nh));
                    for headingIndex = 1:nh
                        [phi_I(:,headingIndex),dphi_dn] = decompose_incident_wave_symmetry(body_centers, ...
                            body_normals,omega,cfg.grav,cfg.water_depth,cfg.wave.headings(headingIndex),1, ...
                            cfg.isx,cfg.isy,parity);
% The assembled source operator uses the inward-body
% normal convention; physical q_D=-q_I therefore
% enters this algebraic system with +S*q_I.
                        rhs(:,headingIndex) = S * dphi_dn;
                    end
                    [phi_D,~] = solve_complex_system(K,rhs,factors);
                    entry.phi_incident_components(:,:,symmetryGroupIndex) = phi_I;
                    entry.phi_diff_components(:,:,symmetryGroupIndex) = phi_D(1:nb,:);
                    entry.phi_incident_body = entry.phi_incident_body + phi_I;
                    entry.phi_diff_body = entry.phi_diff_body + phi_D(1:nb,:);
                end
            end
            pot_cache.entries(k) = entry;
        else
            entry = pot_cache.entries(k);
        end

        if cfg.run.iforce == 1
            if cfg.run.irad == 1
                [A_all(:,:,k),pressure_damping,coefficient_diagnostics] = compute_hydrodynamic_coeffs( ...
                    entry.phi_rad_body,nj,body_areas,omega,cfg.rho,symmetry);
                [B_all(:,:,k),energy_diagnostics] = compute_radiation_damping_energy( ...
                    entry.phi_rad_body,nj,body_centers,body_normals,body_areas, ...
                    omega,cfg,mode_parity,72);
                coefficient_diagnostics.pressure_damping = pressure_damping;
                coefficient_diagnostics.energy_damping = B_all(:,:,k);
                coefficient_diagnostics.energy_diagnostics = energy_diagnostics;
                coefficient_diagnostics.damping_symmetry_error = norm(B_all(:,:,k) - B_all(:,:,k).','fro') / ...
                    max(norm(B_all(:,:,k),'fro'),eps);
                coefficient_diagnostics.min_damping_diagonal = min(diag(B_all(:,:,k)));
                diagnostics(k) = coefficient_diagnostics;
            end
            if cfg.run.idiff == 1
                weights = symmetry_force_weights(mode_parity,group_parity,cfg.isx,cfg.isy);
                for symmetryGroupIndex = 1:ng
                    [component_force,component_pressure] = compute_wave_excitation( ...
                        entry.phi_incident_components(:,:,symmetryGroupIndex), ...
                        entry.phi_diff_components(:,:,symmetryGroupIndex), ...
                        nj,body_areas,omega,cfg.rho);
                    F_all(:,:,k) = F_all(:,:,k) + ...
                        component_force .* weights(:,symmetryGroupIndex);
                    pressure_all(:,:,k) = pressure_all(:,:,k) + component_pressure;
                end
            end
        end
        fprintf('[OK] Frequency %d/%d completed.\n', k, nf);
    end
    if cfg.run.ipoten == 1
        save_potential_cache(cfg.files.potential_cache,pot_cache);
    end

    %% Stage 5: Solve response and second-order mean-drift loads

    rao = [];
    if cfg.run.iforce == 1 && cfg.run.irad == 1 && cfg.run.idiff == 1
        rao = solve_rao(cfg.freq.omegas,A_all,B_all,C,F_all,cfg);
    end
    drift = struct('enabled',false,'near',[],'far',[],'near_details',{{}}, ...
'far_details',{{}},'csv_file','','mat_file','');
    if cfg.run.idrift == 1
        if isempty(rao)
            error('CRESTU:DriftPrerequisites', ...
'IDRIFT=1 requires IFORCE=IRAD=IDIFF=1 and a solved/loaded potential cache.');
        end
        drift = compute_case_drift(domain,pot_cache,rao,F_all,nj, ...
            mode_parity,group_parity);
    end
    %% Stage 6: Preserve result fields and export the result file

    results = struct('schema_version',5,'case_name',cfg.case_name,'config',cfg,'stats',domain.stats, ...
'omegas',cfg.freq.omegas,'headings',cfg.wave.headings,'symmetry',symmetry, ...
'generalized_normals',nj,'added_mass',A_all,'damping',B_all,'excitation',F_all, ...
'pressure',pressure_all,'hydrostatic_matrix',C,'hydrostatics',hydrostatics, ...
'rao',rao,'drift',drift,'diagnostics',diagnostics,'assembly_info',{assembly_info}, ...
'potential_cache_file',cfg.files.potential_cache);
    save(cfg.files.results,'results','-v7.3');
    fprintf('[OK] CRESTU case completed | results = %s\n', cfg.files.results);
end

function frequency_domain = build_frequency_local_domain(base_domain,omega)
% BUILD_FREQUENCY_LOCAL_DOMAIN Build frequency-dependent outer boundaries and merged geometry.
%
% Syntax:
%   frequency_domain=build_frequency_local_domain(base_domain,omega)
%
% Description:
%   Build frequency-dependent outer boundaries and merged geometry.
%   The implementation preserves the CRESTU solver and data-field contracts.
%
% Inputs:
%   base_domain - Frequency-independent boundary-domain structure [-].
%   omega - Wave angular frequency [rad/s].
%
% Outputs:
%   frequency_domain - Frequency-local boundary-domain structure [-].
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
%BUILD_FREQUENCY_LOCAL_DOMAIN Rebuild outer boundaries for one wavelength.
%% Stage 1: Initialize inputs and dependencies

    cfg = base_domain.cfg;
    cfg.fs.r_outer = cfg.fs.tuning_report.original_outer_radius;
    [cfg,cfg.fs.tuning_report] = tune_sponge_layer(cfg,omega);

%% Stage 2: Run the core calculation

    waterlines = cell(cfg.n_bodies,1);
    for body_index = 1:cfg.n_bodies
        waterlines{body_index} = extract_waterline(base_domain.body_list{body_index},cfg.z_tol);
    end
    if cfg.n_bodies == 1
        mesh_fs = generate_free_surface_bmf(cfg.files.fs,waterlines{1},cfg);
    else
        mesh_fs = generate_multibody_free_surface_bmf(cfg.files.fs,waterlines,cfg);
    end
    mesh_seabed = [];
    mesh_farfield = [];
    if cfg.water_depth>0
        reference_waterline = waterlines{1};
        if (cfg.isx || cfg.isy) && cfg.n_bodies == 1
            full_waterline = complete_waterline_by_symmetry(reference_waterline,cfg.isx,cfg.isy);
            cfg_full = cfg;
            cfg_full.isx = 0;
            cfg_full.isy = 0;
            mesh_seabed = generate_reduced_seabed_mesh(cfg,full_waterline.n_pts);
            mesh_farfield = generate_farfield_mesh(full_waterline,cfg_full,cfg.fs.nz_farfield);
            mesh_farfield = reduce_mesh_by_symmetry(mesh_farfield,cfg.isx,cfg.isy,cfg.z_tol);
        elseif cfg.n_bodies == 1
            mesh_seabed = generate_seabed_mesh(reference_waterline,cfg,mesh_fs);
            mesh_farfield = generate_farfield_mesh(reference_waterline,cfg,cfg.fs.nz_farfield);
        else
            mesh_seabed = generate_seabed_disk_mesh(cfg,mesh_fs);
            mesh_farfield = generate_farfield_mesh(reference_waterline,cfg,cfg.fs.nz_farfield);
        end
        write_bmf(cfg.files.seabed,mesh_seabed);
        write_bmf(cfg.files.farfield,mesh_farfield);
    end
    stats = struct('total_body_panels',base_domain.stats.total_body_panels, ...
'fs_panels',mesh_fs.n_panels,'seabed_panels',local_panel_count(mesh_seabed), ...
'farfield_panels',local_panel_count(mesh_farfield));
    stats.total_dofs = stats.total_body_panels + stats.fs_panels + stats.seabed_panels + stats.farfield_panels;
    frequency_domain = struct('cfg',cfg,'body_list',{base_domain.body_list},'fs',mesh_fs, ...
'seabed',mesh_seabed,'farfield',mesh_farfield,'waterline',{waterlines},'stats',stats);
    frequency_domain.geometry = merge_domain_geometry(frequency_domain);
end

function count = local_panel_count(mesh)
% LOCAL_PANEL_COUNT Return the panel count for an optional boundary mesh.
%
% Syntax:
%   count=local_panel_count(mesh)
%
% Description:
%   Return the panel count for an optional boundary mesh.
%   The implementation preserves the CRESTU solver and data-field contracts.
%
% Inputs:
%   mesh - Panel-mesh structure or empty value [-].
%
% Outputs:
%   count - Panel count, nonnegative integer [-].
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
%LOCAL_PANEL_COUNT Return zero for an absent optional boundary mesh.
%% Stage 1: Initialize inputs and dependencies

    if isempty(mesh)
        count = 0;
    else
        count = mesh.n_panels;
    end
end

function value = empty_diagnostics(ndof)
% EMPTY_DIAGNOSTICS Create an initialized hydrodynamic diagnostics structure.
%
% Syntax:
%   value=empty_diagnostics(ndof)
%
% Description:
%   Create an initialized hydrodynamic diagnostics structure.
%   The implementation preserves the CRESTU solver and data-field contracts.
%
% Inputs:
%   ndof - Number of global rigid-body degrees of freedom [-].
%
% Outputs:
%   value - Initialized diagnostics structure [-].
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

    value = struct('radiation_force',complex(zeros(ndof)),'added_mass_symmetry_error',NaN, ...
'raw_added_mass_symmetry_error',NaN,'raw_damping_symmetry_error',NaN, ...
'damping_symmetry_error',NaN,'min_added_mass_diagonal',NaN, ...
'min_damping_diagonal',NaN,'symmetry_weights',ones(ndof), ...
'pressure_damping',zeros(ndof),'energy_damping',zeros(ndof), ...
'energy_diagnostics',struct());
end
function signature = environment_signature(cfg)
% ENVIRONMENT_SIGNATURE Create a numeric signature for cache compatibility checks.
%
% Syntax:
%   signature=environment_signature(cfg)
%
% Description:
%   Create a numeric signature for cache compatibility checks.
%   The implementation preserves the CRESTU solver and data-field contracts.
%
% Inputs:
%   cfg - Validated CRESTU configuration structure [-].
%
% Outputs:
%   signature - Numeric compatibility signature vector [-].
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

    signature = [cfg.water_depth,cfg.grav,cfg.rho,cfg.fs.r_inner,cfg.fs.r_outer, ...
        cfg.fs.mu0,cfg.isx,cfg.isy];
end
function signature = geometry_signature(g)
% GEOMETRY_SIGNATURE Create a numeric geometry signature for cache compatibility checks.
%
% Syntax:
%   signature=geometry_signature(g)
%
% Description:
%   Create a numeric geometry signature for cache compatibility checks.
%   The implementation preserves the CRESTU solver and data-field contracts.
%
% Inputs:
%   g - Merged panel-geometry structure with coordinates in meters [m].
%
% Outputs:
%   signature - Numeric compatibility signature vector [-].
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

    signature = [sum(g.centers(:)),sum(g.centers(:) .^ 2),sum(g.normals(:)), ...
        sum(g.areas(:)),sum(g.vertices(:)),sum(g.vertices(:) .^ 2)];
end

function drift = compute_case_drift(domain,pot_cache,rao,F_all,nj,mode_parity,group_parity)
% COMPUTE_CASE_DRIFT Reconstruct full fields and evaluate near- and far-field mean drift.
%
% Syntax:
%   drift=compute_case_drift(domain,pot_cache,rao,F_all,nj,mode_parity,group_parity)
%
% Description:
%   Reconstruct full fields and evaluate near- and far-field mean drift.
%   The implementation preserves the CRESTU solver and data-field contracts.
%
% Inputs:
%   domain - Prepared boundary-domain structure [-].
%   pot_cache - Potential-cache structure [-].
%   rao - Complex rigid-body response amplitudes, [m] or [rad].
%   F_all - Complex first-order generalized forces, [N] or [N*m].
%   nj - Generalized body-normal matrix [-] or [m].
%   mode_parity - Parity flags for global rigid-body modes [-].
%   group_parity - Unique parity groups for symmetry solves [-].
%
% Outputs:
%   drift - Near- and far-field mean-drift result structure [-].
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
%COMPUTE_CASE_DRIFT Reconstruct full fields and evaluate near/far mean loads.
% Reconstruct symmetry-reduced solutions by parity before second-order integration.
%% Stage 1: Initialize inputs and dependencies

    cfg = domain.cfg;

%% Stage 2: Run the core calculation

    nf = numel(cfg.freq.omegas);
    nh = numel(cfg.wave.headings);
    ndof = 6 * cfg.n_bodies;
    ng = size(group_parity,1);
    near_loads = zeros(ndof,nh,nf);
    far_loads = zeros(6,nh,nf);
    near_details = cell(nf,cfg.n_bodies);
    far_details = cell(nf,1);
    body_start = 1;
    body_ranges = cell(cfg.n_bodies,1);
    for bodyIndex = 1:cfg.n_bodies
        body_ranges{bodyIndex} = body_start:(body_start + domain.body_list{bodyIndex}.n_panels - 1);
        body_start = body_ranges{bodyIndex}(end) + 1;
    end
    for k = 1:nf
        omega = cfg.freq.omegas(k);
        xi = rao.complex(:,:,k);
        entry = pot_cache.entries(k);
        ensemble_mesh = [];
        ensemble_phi = [];
        ensemble_q = [];
        for bodyIndex = 1:cfg.n_bodies
            reduced = domain.body_list{bodyIndex};
            panel_rows = body_ranges{bodyIndex};
            force_rows = (bodyIndex - 1) * 6 + (1:6);
            full_mesh = expand_mesh_by_symmetry(reduced,cfg.isx,cfg.isy);
            phi_incident = complex(zeros(full_mesh.n_panels,nh));
            phi_diffraction = complex(zeros(full_mesh.n_panels,nh));
            for normalVelocity = 1:ng
                incident_component = entry.phi_incident_components(panel_rows,:,normalVelocity);
                diffraction_component = entry.phi_diff_components(panel_rows,:,normalVelocity);
                phi_incident = phi_incident + expand_scalar_by_symmetry(incident_component, ...
                    repmat(group_parity(normalVelocity,:),nh,1),cfg.isx,cfg.isy);
                phi_diffraction = phi_diffraction + expand_scalar_by_symmetry(diffraction_component, ...
                    repmat(group_parity(normalVelocity,:),nh,1),cfg.isx,cfg.isy);
            end
            phi_modes = expand_scalar_by_symmetry(entry.phi_rad_body(panel_rows,:), ...
                mode_parity,cfg.isx,cfg.isy);
            q_modes = expand_scalar_by_symmetry(1i * omega * nj(panel_rows,:), ...
                mode_parity,cfg.isx,cfg.isy);
            phi_radiation = phi_modes * xi;
            q_radiation = q_modes * xi;
            phi_total = phi_incident + phi_diffraction + phi_radiation;
            q_total = q_radiation;
% Pinkster Term 3 uses excitation only. Radiation is already
% represented through the total first-order potential.
            first_order_force = F_all(force_rows,:,k);
            state = struct('phi',phi_total,'dphi_dn',q_total, ...
'rao',xi(force_rows,:),'first_order_force',first_order_force, ...
'omega',omega,'headings',cfg.wave.headings);
            cfg_full = cfg;
            cfg_full.isx = 0;
            cfg_full.isy = 0;
            waterline = extract_waterline(full_mesh,cfg.z_tol);
            detail = compute_drift_nearfield(full_mesh,waterline,state,cfg_full);
            near_details{k,bodyIndex} = detail;
            near_loads(force_rows,:,k) = detail.total;
            q_incident = complex(zeros(full_mesh.n_panels,nh));
            for heading_index = 1:nh
                [~,q_incident(:,heading_index)] = compute_incident_wave( ...
                    full_mesh.centers,full_mesh.normals,omega,cfg.grav,cfg.water_depth, ...
                    cfg.wave.headings(heading_index),1);
            end
% The Kochin function describes the outgoing disturbance field.
% The incident field is represented by the explicit cross term
% in the Maruo-Newman force formula and must not be integrated
% a second time as part of H(theta).
            phi_disturbance = phi_diffraction + phi_radiation;
            q_disturbance = -q_incident + q_radiation;
            [ensemble_mesh,ensemble_phi,ensemble_q] = append_body_state( ...
                ensemble_mesh,ensemble_phi,ensemble_q,full_mesh,phi_disturbance,q_disturbance);
        end
        far_state = struct('phi',ensemble_phi,'dphi_dn',ensemble_q,'omega',omega, ...
'headings',cfg.wave.headings);
        cfg_full = cfg;
        cfg_full.isx = 0;
        cfg_full.isy = 0;
        far_details{k} = compute_drift_farfield(ensemble_mesh,far_state,cfg_full);
        far_loads(:,:,k) = far_details{k}.total;
    end
    characteristic_length = 0;
    for bodyIndex = 1:cfg.n_bodies
        body_vertices = expand_mesh_by_symmetry(domain.body_list{bodyIndex},cfg.isx,cfg.isy).vertices;
        body_extent = zeros(1,3);
        for coordinate = 1:3
            coordinate_values = reshape(body_vertices(:,:,coordinate),[],1);
            body_extent(coordinate) = max(coordinate_values) - min(coordinate_values);
        end
        characteristic_length = max(characteristic_length,max(body_extent));
    end
    csv_file = fullfile(cfg.config_dir,sprintf('%s_MeanDrift.csv',cfg.case_name));
    mat_file = fullfile(cfg.config_dir,sprintf('%s_MeanDrift.mat',cfg.case_name));
    export_drift_loads(csv_file,cfg.freq.omegas,cfg.wave.headings,near_loads,far_loads, ...
        cfg.rho,cfg.grav,1,characteristic_length);
    export_drift_loads(mat_file,cfg.freq.omegas,cfg.wave.headings,near_loads,far_loads, ...
        cfg.rho,cfg.grav,1,characteristic_length);
    drift = struct('enabled',true,'near',near_loads,'far',far_loads, ...
'near_details',{near_details},'far_details',{far_details}, ...
'normalization',0.5 * cfg.rho * cfg.grav * characteristic_length, ...
'wave_amplitude',1,'characteristic_length',characteristic_length, ...
'csv_file',csv_file,'mat_file',mat_file);
end

function [combinedMesh, velocityPotential, normalVelocity] = ...
    append_body_state(combinedMesh, velocityPotential, normalVelocity, ...
    bodyMesh, newVelocityPotential, newNormalVelocity)
% APPEND_BODY_STATE Append one body mesh and its potential state to an ensemble.
%
% Syntax:
%   [combinedMesh, velocityPotential, normalVelocity] = ...
%       append_body_state(combinedMesh, velocityPotential, normalVelocity, ...
%       bodyMesh, newVelocityPotential, newNormalVelocity)
%
% Description:
%   Append one body mesh and its potential state to an ensemble.
%   The implementation preserves the CRESTU solver and data-field contracts.
%
% Inputs:
%   combined - Combined body-mesh structure or empty value [-].
%   phi - Combined velocity-potential array [m^2/s].
%   normalVelocity - Combined normal velocity array [m/s].
%   mesh - Panel-mesh structure or empty value [-].
%   new_phi - Velocity potential for the appended body [m^2/s].
%   new_q - Normal velocity for the appended body [m/s].
%
% Outputs:
%   combinedMesh - Updated combined body-mesh structure [-].
%   velocityPotential - Updated combined velocity-potential array [m^2/s].
%   normalVelocity - Updated combined normal-velocity array [m/s].
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

    if isempty(combinedMesh)
        combinedMesh = bodyMesh;
        velocityPotential = newVelocityPotential;
        normalVelocity = newNormalVelocity;
        return;
    end

%% Stage 2: Append body arrays to the ensemble

    fields = {'vertices','centers','normals','areas','e1','e2','panel_type'};
    for fieldIndex = 1:numel(fields)
        name = fields{fieldIndex};
        combinedMesh.(name) = cat(1, combinedMesh.(name), bodyMesh.(name));
    end
    combinedMesh.n_panels = size(combinedMesh.centers, 1);
    velocityPotential = [velocityPotential; newVelocityPotential];
    normalVelocity = [normalVelocity; newNormalVelocity];
end
