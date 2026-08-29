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

    %% 阶段 3: 建立完整缓存清单并初始化势函数缓存

    baseMeshAudit = audit_domain_meshes(domain);
    codeVersion = get_rankine_code_version();
    [baseCacheKey, baseCacheSpecification] = build_base_cache_key( ...
        cfg, baseMeshAudit, mode_parity, group_parity, codeVersion);
    emptyEntry = create_empty_cache_entry(nb, ndof, nh, ng);

    if cfg.run.ipoten == 0
        pot_cache = load_potential_cache(cfg.files.potential_cache, ...
            baseCacheKey);
    else
        fprintf('[CACHE] MISS | case=%s | reason=forced_recompute\n', ...
            cfg.case_name);
        pot_cache = struct('schema_version', 5, ...
            'case_name', cfg.case_name, 'omegas', cfg.freq.omegas, ...
            'headings', cfg.wave.headings, 'n_body_panels', nb, ...
            'n_total_panels', geom.total_panels, ...
            'mode_parity', mode_parity, 'group_parity', group_parity, ...
            'n_dofs', ndof, 'base_cache_key', baseCacheKey, ...
            'base_cache_specification', baseCacheSpecification, ...
            'code_version', codeVersion, ...
            'created', char(datetime('now', ...
            'Format', 'yyyyMMdd''T''HHmmss')), ...
            'entries', repmat(emptyEntry, nf, 1));
    end

    A_all = zeros(ndof,ndof,nf);
    B_all = zeros(ndof,ndof,nf);
    A_raw_all = zeros(ndof,ndof,nf);
    B_raw_all = zeros(ndof,ndof,nf);
    B_energy_all = zeros(ndof,ndof,nf);
    F_all = complex(zeros(ndof,nh,nf));
    pressure_all = complex(zeros(nb,nh,nf));
    diagnostics = repmat(empty_diagnostics(ndof),nf,1);
    assembly_info = cell(nf,ng);
    auditEntries = repmat(empty_audit_entry(), nf, 1);

    %% 阶段 4: 按频率重建真实边界、验证缓存并求解一阶势

    for k = 1:nf
        omega = cfg.freq.omegas(k); % [rad/s]
        fprintf('[INFO] Frequency %d/%d | omega = %.6g rad/s\n', k, nf, omega);
        frequency_domain = build_frequency_local_domain(domain,omega);
        frequency_geometry = frequency_domain.geometry;
        frequency_cfg = frequency_domain.cfg;
        frequencyMeshAudit = audit_domain_meshes(frequency_domain);
        [entryCacheKey, entryCacheSpecification] = build_potential_cache_key( ...
            frequency_cfg, frequencyMeshAudit, omega, mode_parity, ...
            group_parity, codeVersion);

        if cfg.run.ipoten == 1
            entry = pot_cache.entries(k);
            entry.omega = omega;
            entry.frequency_domain_stats = frequency_domain.stats;
            entry.frequency_environment = [frequency_cfg.fs.r_inner,frequency_cfg.fs.r_outer, ...
                frequency_cfg.fs.tuning_report.dynamic_mu0];
            entry.cache_key = entryCacheKey;
            entry.cache_specification = entryCacheSpecification;
            entry.mesh_audit = frequencyMeshAudit;
            entry.assembly_info = cell(ng, 1);
            entry.solver_diagnostics = cell(ng, 1);
            cacheStatus = 'miss';

            for symmetryGroupIndex = 1:ng
                parity = group_parity(symmetryGroupIndex, :);
                modes = group_modes{symmetryGroupIndex};
                [K,S,assembly_info{k,symmetryGroupIndex}] = assemble_rankine_matrix(frequency_geometry.total_panels, ...
                    frequency_geometry.centers,frequency_geometry.normals,frequency_geometry.vertices, ...
                    frequency_domain.stats,omega,frequency_cfg.grav,frequency_cfg,parity);
                entry.assembly_info{symmetryGroupIndex} = ...
                    assembly_info{k,symmetryGroupIndex};
                factors = [];
                radiationResidual = NaN;
                diffractionResidual = NaN;

                if cfg.run.irad == 1
                    rhs = S * (1i * omega * nj(:,modes));
                    [phi,factors] = solve_complex_system(K,rhs,factors);
                    entry.phi_rad_body(:,modes) = phi(1:nb,:);
                    radiationResidual = relative_linear_residual(K, phi, rhs);
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
                    diffractionResidual = relative_linear_residual(K, phi_D, rhs);
                    entry.phi_incident_components(:,:,symmetryGroupIndex) = phi_I;
                    entry.phi_diff_components(:,:,symmetryGroupIndex) = phi_D(1:nb,:);
                    entry.phi_incident_body = entry.phi_incident_body + phi_I;
                    entry.phi_diff_body = entry.phi_diff_body + phi_D(1:nb,:);
                end

                entry.solver_diagnostics{symmetryGroupIndex} = struct( ...
                    'radiationRelativeResidual', radiationResidual, ...
                    'diffractionRelativeResidual', diffractionResidual);
            end
            pot_cache.entries(k) = entry;
        else
            entry = pot_cache.entries(k);
            validate_cache_entry(entry, entryCacheKey, omega, ...
                frequencyMeshAudit);
            cacheStatus = 'hit';

            for symmetryGroupIndex = 1:ng
                assembly_info{k,symmetryGroupIndex} = ...
                    entry.assembly_info{symmetryGroupIndex};
            end
        end

        auditEntries(k) = build_frequency_audit_entry(cfg.case_name, ...
            omega, frequencyMeshAudit, entryCacheKey, cacheStatus, entry);
        print_frequency_audit(auditEntries(k));

        if cfg.run.iforce == 1
            if cfg.run.irad == 1
                [A_all(:,:,k),pressure_damping,coefficient_diagnostics] = compute_hydrodynamic_coeffs( ...
                    entry.phi_rad_body,nj,body_areas,omega,cfg.rho,symmetry);
                B_all(:,:,k) = pressure_damping;
                A_raw_all(:,:,k) = coefficient_diagnostics.added_mass_raw;
                B_raw_all(:,:,k) = coefficient_diagnostics.damping_raw;
                [B_energy_all(:,:,k),energy_diagnostics] = compute_radiation_damping_energy( ...
                    entry.phi_rad_body,nj,body_centers,body_normals,body_areas, ...
                    omega,cfg,mode_parity,72);
                coefficient_diagnostics.pressure_damping = pressure_damping;
                coefficient_diagnostics.energy_damping = B_energy_all(:,:,k);
                coefficient_diagnostics.energy_diagnostics = energy_diagnostics;
                coefficient_diagnostics.pressure_energy_relative_residual = ...
                    relative_matrix_difference(B_all(:,:,k), B_energy_all(:,:,k));
                coefficient_diagnostics.B33_pressure_energy_relative_residual = ...
                    abs(B_all(3,3,k) - B_energy_all(3,3,k)) / ...
                    max([abs(B_all(3,3,k)), abs(B_energy_all(3,3,k)), eps]);
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

    %% 阶段 5: 求解刚体响应并按开关处理二阶平均漂移

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
    %% 阶段 6: 保存原始矩阵、互易投影、能量核查与审计清单

    meshAuditFile = export_mesh_audit(cfg, auditEntries);
    auditSummary = struct('baseMesh', baseMeshAudit, ...
        'frequencyEntries', auditEntries, ...
        'baseCacheKey', baseCacheKey, ...
        'baseCacheSpecification', baseCacheSpecification, ...
        'codeVersion', codeVersion, ...
        'meshAuditFile', meshAuditFile, ...
        'cacheMode', cache_mode_name(cfg.run.ipoten));
    results = struct('schema_version',6,'case_name',cfg.case_name,'config',cfg,'stats',domain.stats, ...
'omegas',cfg.freq.omegas,'headings',cfg.wave.headings,'symmetry',symmetry, ...
'generalized_normals',nj,'A_raw',A_raw_all,'B_raw',B_raw_all, ...
'A_reciprocal',A_all,'B_reciprocal',B_all, ...
'added_mass_raw',A_raw_all,'damping_raw',B_raw_all, ...
'added_mass',A_all,'damping',B_all,'damping_energy',B_energy_all, ...
'excitation',F_all, ...
'pressure',pressure_all,'hydrostatic_matrix',C,'hydrostatics',hydrostatics, ...
'rao',rao,'drift',drift,'diagnostics',diagnostics,'assembly_info',{assembly_info}, ...
'potential_cache_file',cfg.files.potential_cache,'audit',auditSummary);
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

    bodyWaterlines = base_domain.body_waterlines;
    outerWaterlines = base_domain.outer_waterlines;

    if cfg.n_bodies == 1
        mesh_fs = generate_free_surface_bmf(cfg.files.fs, ...
            outerWaterlines{1}, cfg);
    else
        mesh_fs = generate_multibody_free_surface_bmf(cfg.files.fs, ...
            outerWaterlines, cfg);
    end
    mesh_seabed = [];
    mesh_farfield = [];
    if cfg.water_depth>0
        reference_waterline = outerWaterlines{1};
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
'body_vertices',base_domain.stats.body_vertices, ...
'waterline_panels',base_domain.stats.waterline_panels, ...
'outer_waterline_panels',base_domain.stats.outer_waterline_panels, ...
'fs_panels',mesh_fs.n_panels,'seabed_panels',local_panel_count(mesh_seabed), ...
'farfield_panels',local_panel_count(mesh_farfield));
    stats.total_dofs = stats.total_body_panels + stats.fs_panels + stats.seabed_panels + stats.farfield_panels;
    frequency_domain = struct('cfg',cfg,'body_list',{base_domain.body_list},'fs',mesh_fs, ...
'seabed',mesh_seabed,'farfield',mesh_farfield, ...
'waterline',{bodyWaterlines},'body_waterlines',{bodyWaterlines}, ...
'outer_waterlines',{outerWaterlines},'stats',stats);
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

    value = struct('radiation_force',complex(zeros(ndof)), ...
'added_mass_raw',zeros(ndof),'damping_raw',zeros(ndof), ...
'added_mass_reciprocal',zeros(ndof),'damping_reciprocal',zeros(ndof), ...
'raw_added_mass_symmetry_error',NaN,'raw_damping_symmetry_error',NaN, ...
'added_mass_symmetry_error',NaN,'damping_symmetry_error',NaN, ...
'min_added_mass_diagonal',NaN,'min_damping_diagonal',NaN, ...
'symmetry_weights',ones(ndof),'pressure_damping',zeros(ndof), ...
'energy_damping',zeros(ndof),'energy_diagnostics',struct(), ...
'pressure_energy_relative_residual',NaN, ...
'B33_pressure_energy_relative_residual',NaN);
end

function entry = create_empty_cache_entry(bodyPanelCount, globalDofCount, ...
        headingCount, symmetryGroupCount)
% CREATE_EMPTY_CACHE_ENTRY Initialize one schema-5 frequency cache entry.

    entry = struct('omega', 0.0, ... % [rad/s]
        'phi_rad_body', complex(zeros(bodyPanelCount, globalDofCount)), ...
        'phi_diff_body', complex(zeros(bodyPanelCount, headingCount)), ...
        'phi_incident_body', complex(zeros(bodyPanelCount, headingCount)), ...
        'phi_diff_components', complex(zeros(bodyPanelCount, ...
        headingCount, symmetryGroupCount)), ...
        'phi_incident_components', complex(zeros(bodyPanelCount, ...
        headingCount, symmetryGroupCount)), ...
        'frequency_domain_stats', struct(), ...
        'frequency_environment', zeros(1, 3), ...
        'cache_key', '', 'cache_specification', struct(), ...
        'mesh_audit', struct(), ...
        'assembly_info', {cell(symmetryGroupCount, 1)}, ...
        'solver_diagnostics', {cell(symmetryGroupCount, 1)});
end

function [cacheKey, cacheSpecification] = build_base_cache_key( ...
        cfg, meshAudit, modeParity, groupParity, codeVersion)
% BUILD_BASE_CACHE_KEY Fingerprint the complete frequency-independent case state.

    %% 阶段 1: 整理质量、惯量与求解选项

    massSignature = zeros(cfg.n_bodies, 13);

    for bodyIndex = 1:cfg.n_bodies
        massSignature(bodyIndex, :) = [cfg.mass_props(bodyIndex).mass, ...
            cfg.mass_props(bodyIndex).cg, ...
            reshape(cfg.mass_props(bodyIndex).inertia, 1, 9)];
    end

    sourceOrientation = get_rankine_source_orientation();
    solverOptions = struct('potentialProblemFlags', ...
        [cfg.run.irad, cfg.run.idiff], ...
        'waterDepth', cfg.water_depth, 'gravity', cfg.grav, ...
        'density', cfg.rho, 'symmetryFlags', [cfg.isx, cfg.isy], ...
        'waterlineTolerance', cfg.z_tol, ...
        'freeSurfaceParameters', [cfg.fs.nr_near, cfg.fs.nr_sponge, ...
        cfg.fs.r_inner, cfg.fs.r_outer, cfg.fs.sponge_ratio, ...
        cfg.fs.nz_farfield, cfg.fs.mu0], ...
        'doubleLayerSourceNormalConvention', ...
        char(sourceOrientation.convention), ...
        'doubleLayerComponentSigns', sourceOrientation.componentSigns, ...
        'doubleLayerOrientationHash', sourceOrientation.signatureHash, ...
        'massProperties', massSignature);

    %% 阶段 2: 构造并哈希基础缓存清单

    cacheSpecification = struct('schemaVersion', 1, ...
        'caseName', cfg.case_name, ...
        'bodyMeshHash', meshAudit.body.hash, ...
        'freeSurfaceMeshHash', meshAudit.freeSurface.hash, ...
        'bottomMeshHash', meshAudit.bottom.hash, ...
        'outerBoundaryMeshHash', meshAudit.outerBoundary.hash, ...
        'frequencyGrid', cfg.freq.omegas, ... % [rad/s]
        'headings', cfg.wave.headings, ... % [deg]
        'globalDofCount', 6 * cfg.n_bodies, ...
        'calculatedLocalDofs', cfg.calc_modes, ...
        'problemType', sprintf('radiation=%d;diffraction=%d', ...
        cfg.run.irad, cfg.run.idiff), ...
        'modeParity', modeParity, 'groupParity', groupParity, ...
        'solverOptions', solverOptions, ...
        'codeVersion', codeVersion.fingerprint);
    cacheKey = sha256_hash(jsonencode(cacheSpecification));
end

function validate_cache_entry(entry, expectedCacheKey, omega, meshAudit)
% VALIDATE_CACHE_ENTRY Require exact per-frequency key and matrix metadata.

    %% 阶段 1: 检查逐频率缓存字段

    requiredFields = {'omega', 'cache_key', 'mesh_audit', ...
        'assembly_info', 'solver_diagnostics'};

    for fieldIndex = 1:numel(requiredFields)
        if ~isfield(entry, requiredFields{fieldIndex})
            fprintf('[CACHE] MISS | omega=%.9g | reason=entry_field_%s\n', ...
                omega, requiredFields{fieldIndex});
            error('CRESTU:CacheSchema', ...
                'Frequency cache entry lacks %s.', requiredFields{fieldIndex});
        end
    end

    %% 阶段 2: 对比频率、网格和完整缓存键

    if abs(entry.omega - omega) > 8 * eps(max(1, abs(omega)))
        fprintf('[CACHE] MISS | omega=%.9g | reason=frequency_mismatch\n', ...
            omega);
        error('CRESTU:CacheMismatch', ...
            'Frequency cache entry does not match omega=%.17g rad/s.', omega);
    end

    if ~strcmp(entry.cache_key, expectedCacheKey)
        fprintf('[CACHE] MISS | omega=%.9g | reason=entry_key_mismatch | expected=%s\n', ...
            omega, expectedCacheKey);
        error('CRESTU:CacheMismatch', ...
            'Per-frequency potential cache key mismatch.');
    end

    if ~strcmp(entry.mesh_audit.geometry.hash, meshAudit.geometry.hash)
        fprintf('[CACHE] MISS | omega=%.9g | reason=geometry_hash_mismatch\n', ...
            omega);
        error('CRESTU:CacheMismatch', ...
            'Cached merged-geometry hash does not match the active mesh.');
    end

    for groupIndex = 1:numel(entry.assembly_info)
        assemblyRecord = entry.assembly_info{groupIndex};

        if ~isstruct(assemblyRecord) || ...
                ~isfield(assemblyRecord, 'n_unknowns') || ...
                assemblyRecord.n_unknowns ~= meshAudit.bemUnknownCount
            fprintf('[CACHE] MISS | omega=%.9g | reason=matrix_size_mismatch\n', ...
                omega);
            error('CRESTU:CacheMismatch', ...
                'Cached influence-matrix dimension does not match the active mesh.');
        end
    end

    fprintf('[CACHE] HIT | omega=%.9g | key=%s\n', omega, expectedCacheKey);
end

function residual = relative_linear_residual(systemMatrix, solution, rightHandSide)
% RELATIVE_LINEAR_RESIDUAL Evaluate ||K*phi-rhs||_F / ||rhs||_F.

    residual = norm(systemMatrix * solution - rightHandSide, 'fro') / ...
        max(norm(rightHandSide, 'fro'), eps);
end

function residual = relative_matrix_difference(firstMatrix, secondMatrix)
% RELATIVE_MATRIX_DIFFERENCE Compare two matrices without hiding either result.

    residual = norm(firstMatrix - secondMatrix, 'fro') / ...
        max([norm(firstMatrix, 'fro'), norm(secondMatrix, 'fro'), eps]);
end

function auditEntry = empty_audit_entry()
% EMPTY_AUDIT_ENTRY Initialize a scalar record for CSV-safe assignment.

    auditEntry = struct('caseName', '', 'omegaRadPerSecond', 0.0, ...
        'bodyPanelCount', 0, 'bodyVertexCount', 0, ...
        'waterlinePanelCount', 0, 'freeSurfacePanelCount', 0, ...
        'bottomPanelCount', 0, 'outerBoundaryPanelCount', 0, ...
        'bemUnknownCount', 0, 'influenceMatrixRows', 0, ...
        'influenceMatrixColumns', 0, 'sourcePanelCount', 0, ...
        'collocationPointCount', 0, 'bodyCharacteristicSizeM', NaN, ...
        'outerCharacteristicSizeM', NaN, ...
        'maxRadiationLinearResidual', NaN, ...
        'maxDiffractionLinearResidual', NaN, ...
        'bodyMeshHash', '', 'freeSurfaceMeshHash', '', ...
        'bottomMeshHash', '', 'outerBoundaryMeshHash', '', ...
        'bodyFreeSurfaceMeshHash', '', 'sourcePointHash', '', ...
        'collocationPointHash', '', 'cacheKey', '', 'cacheStatus', '');
end

function auditEntry = build_frequency_audit_entry(caseName, omega, ...
        meshAudit, cacheKey, cacheStatus, cacheEntry)
% BUILD_FREQUENCY_AUDIT_ENTRY Capture actual meshes and matrix dimensions.

    unknownCount = meshAudit.bemUnknownCount;
    outerPanelCount = meshAudit.freeSurface.panelCount + ...
        meshAudit.bottom.panelCount + meshAudit.outerBoundary.panelCount;
    outerArea = meshAudit.freeSurface.surfaceArea + ... % [m^2]
        meshAudit.bottom.surfaceArea + ...
        meshAudit.outerBoundary.surfaceArea; % [m^2]
    outerCharacteristicSize = sqrt(outerArea / outerPanelCount); % [m]
    radiationResiduals = cellfun(@(item) ...
        item.radiationRelativeResidual, cacheEntry.solver_diagnostics);
    diffractionResiduals = cellfun(@(item) ...
        item.diffractionRelativeResidual, cacheEntry.solver_diagnostics);
    auditEntry = struct('caseName', caseName, ...
        'omegaRadPerSecond', omega, ... % [rad/s]
        'bodyPanelCount', meshAudit.body.panelCount, ...
        'bodyVertexCount', meshAudit.body.vertexCount, ...
        'waterlinePanelCount', meshAudit.bodyWaterlinePanelCount, ...
        'freeSurfacePanelCount', meshAudit.freeSurface.panelCount, ...
        'bottomPanelCount', meshAudit.bottom.panelCount, ...
        'outerBoundaryPanelCount', meshAudit.outerBoundary.panelCount, ...
        'bemUnknownCount', unknownCount, ...
        'influenceMatrixRows', unknownCount, ...
        'influenceMatrixColumns', unknownCount, ...
        'sourcePanelCount', meshAudit.sourcePanelCount, ...
        'collocationPointCount', meshAudit.collocationPointCount, ...
        'bodyCharacteristicSizeM', meshAudit.body.characteristicSize, ... % [m]
        'outerCharacteristicSizeM', outerCharacteristicSize, ... % [m]
        'maxRadiationLinearResidual', max(radiationResiduals), ...
        'maxDiffractionLinearResidual', max(diffractionResiduals), ...
        'bodyMeshHash', meshAudit.body.hash, ...
        'freeSurfaceMeshHash', meshAudit.freeSurface.hash, ...
        'bottomMeshHash', meshAudit.bottom.hash, ...
        'outerBoundaryMeshHash', meshAudit.outerBoundary.hash, ...
        'bodyFreeSurfaceMeshHash', meshAudit.bodyFreeSurfaceHash, ...
        'sourcePointHash', meshAudit.sourcePointHash, ...
        'collocationPointHash', meshAudit.collocationPointHash, ...
        'cacheKey', cacheKey, 'cacheStatus', cacheStatus);
end

function print_frequency_audit(auditEntry)
% PRINT_FREQUENCY_AUDIT Print a structured actual-mesh audit line.

    fprintf(['[AUDIT] case=%s omega=%.6g | body=%d panels/%d vertices ', ...
        '| waterline=%d | fs=%d bottom=%d outer=%d | unknowns=%d ', ...
        '| matrix=%dx%d | cache=%s\n'], ...
        auditEntry.caseName, auditEntry.omegaRadPerSecond, ...
        auditEntry.bodyPanelCount, auditEntry.bodyVertexCount, ...
        auditEntry.waterlinePanelCount, auditEntry.freeSurfacePanelCount, ...
        auditEntry.bottomPanelCount, auditEntry.outerBoundaryPanelCount, ...
        auditEntry.bemUnknownCount, auditEntry.influenceMatrixRows, ...
        auditEntry.influenceMatrixColumns, upper(auditEntry.cacheStatus));
    fprintf('[AUDIT] body_hash=%s | fs_hash=%s | cache_key=%s\n', ...
        auditEntry.bodyMeshHash, auditEntry.freeSurfaceMeshHash, ...
        auditEntry.cacheKey);
end

function auditFilename = export_mesh_audit(cfg, auditEntries)
% EXPORT_MESH_AUDIT Save the per-frequency forensic audit as CSV.

    auditFilename = fullfile(cfg.config_dir, ...
        sprintf('%s_MeshAudit.csv', cfg.case_name));
    auditTable = struct2table(auditEntries);
    writetable(auditTable, auditFilename);
    fprintf('[OK] Mesh/cache audit saved: %s\n', auditFilename);
end

function modeName = cache_mode_name(ipotenFlag)
% CACHE_MODE_NAME Convert the potential control flag to an audit label.

    if ipotenFlag == 1
        modeName = 'clean-recompute';
    else
        modeName = 'cache-reload';
    end
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
