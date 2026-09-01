function summary = Run_Phase3_2B_Legacy_Regression(suitePosition)
% RUN_PHASE3_2B_LEGACY_REGRESSION Run start/end legacy and mode controls.

    arguments
        suitePosition (1, 1) string {mustBeMember(suitePosition, ...
            ["START", "END"])}
    end
    definition = phase32b.initialize();
    outputFile = fullfile(definition.artifactDirectory, ...
        'Phase3_2B_Legacy_Regression.csv');
    rows = repmat(empty_row(), 0, 1);

    %% Stage 1: Rebuild five q-only LEGACY geometry references

    referenceFile = fullfile(definition.convergenceDirectory, ...
        'Phase3_1_Artifacts', 'Phase3_1_Validation_Result.mat');
    loadedReference = load(referenceFile, 'summary');
    reference = loadedReference.summary.radiusTable;
    reference = reference(reference.radius_factor == 1.5, :);
    assert(height(reference) == 5, 'CRESTU:Phase32BLegacyReference', ...
        'Five frozen q-only geometry rows are required.');
    baseline = read_config(definition.baselineConfig);
    formal = read_config(definition.formalConfig);
    for referenceIndex = 1:height(reference)
        omega = reference.omega_rad_s(referenceIndex); % [rad/s]
        runtimeDirectory = tempname(definition.artifactDirectory);
        mkdir(runtimeDirectory);
        cleanup = onCleanup(@() phase32.remove_runtime( ...
            definition, runtimeDirectory));
        caseName = sprintf('P32B_LegacyGeometry_w%04d', round(1000 * omega));
        configFile = fullfile(runtimeDirectory, [caseName, '.cfg']);
        write_q_only_config(configFile, caseName, omega, ...
            reference.q_top(referenceIndex), ...
            reference.q_bottom(referenceIndex), baseline, formal);
        config = read_config(configFile);
        assert(get_polar_topology_mode(config) == "LEGACY", ...
            'CRESTU:Phase32BLegacyMode', ...
            'A q-only case without mode must remain LEGACY.');
        domain = build_bmf_domain(configFile);
        audit = audit_domain_meshes(domain);
        exact = string(audit.body.hash) == ...
            reference.body_geometry_hash(referenceIndex) && ...
            string(audit.freeSurface.hash) == ...
            reference.free_surface_geometry_hash(referenceIndex) && ...
            string(audit.bottom.hash) == ...
            reference.bottom_geometry_hash(referenceIndex) && ...
            string(audit.outerBoundary.hash) == ...
            reference.farfield_geometry_hash(referenceIndex) && ...
            string(audit.geometry.hash) == ...
            reference.combined_geometry_hash(referenceIndex);
        row = empty_row();
        row.suite_position = suitePosition;
        row.test_id = "LEGACY_Q_ONLY_GEOMETRY";
        row.omega_rad_s = omega;
        row.topology_request = "OMITTED";
        row.observed_topology_mode = "LEGACY";
        row.expected_combined_hash = ...
            reference.combined_geometry_hash(referenceIndex);
        row.actual_combined_hash = string(audit.geometry.hash);
        row.exact_hash_match = exact;
        row.status = pass_fail(exact);
        row.notes = "fresh q-only mesh; no explicit topology mode";
        rows(end + 1, 1) = row; %#ok<AGROW>
        clear cleanup
    end

    %% Stage 2: Run default-negative and explicit-V2 positive controls

    negative = empty_row();
    negative.suite_position = suitePosition;
    negative.test_id = "TOPOLOGY_DEFAULT_NEGATIVE_CONTROL";
    negative.topology_request = "OMITTED";
    negative.observed_topology_mode = get_polar_topology_mode(struct());
    negative.status = pass_fail(negative.observed_topology_mode == "LEGACY");
    negative.notes = "missing mode must never select V2";
    rows(end + 1, 1) = negative;

    positive = run_positive_control(definition, suitePosition);
    rows(end + 1, 1) = positive;

    %% Stage 3: Rerun the clean LEGACY 1.5 rad/s response regression

    response = run_legacy_response(definition, suitePosition);
    rows(end + 1, 1) = response;
    current = struct2table(rows);
    if suitePosition == "END" && isfile(outputFile)
        previous = readtable(outputFile, 'TextType', 'string');
        current = [previous; current];
    end
    writetable(current, outputFile);
    selected = current.suite_position == suitePosition;
    pass = nnz(selected) == 8 && all(current.status(selected) == "PASS");
    summary = struct('schemaVersion', 1, ...
        'suitePosition', char(suitePosition), ...
        'testCount', nnz(selected), ...
        'passCount', nnz(current.status(selected) == "PASS"), ...
        'overallPass', pass);
    fprintf('[RESULT] PHASE3_2B_%s_LEGACY_REGRESSION = %s\n', ...
        suitePosition, pass_fail(pass));
    assert(pass, 'CRESTU:Phase32BLegacyRegression', ...
        '%s legacy/control regression did not pass 8/8.', suitePosition);
end

function row = run_positive_control(definition, suitePosition)
    level = definition.levels(3);
    omega = 1.5; % [rad/s]
    expected = phase32b.find_frozen_row(definition, level, omega);
    specification = phase32a.create_case(definition, level, omega, ...
        'phase3_2b_positive_control');
    cleanup = onCleanup(@() phase32.remove_runtime( ...
        definition, specification.runtimeDirectory));
    config = read_config(specification.configFile);
    domain = build_bmf_domain(specification.configFile);
    [actual, ~] = phase32a.audit_case( ...
        definition, specification, domain, NaN);
    [operator, controls, meshAudit] = ...
        phase32b.preview_operator_state(domain, omega);
    phase32b.assert_v5_match(expected, actual, specification, ...
        operator, controls, meshAudit);
    row = empty_row();
    row.suite_position = suitePosition;
    row.test_id = "QUALITY_CONTROLLED_V2_POSITIVE_CONTROL";
    row.omega_rad_s = omega;
    row.topology_request = "EXPLICIT_QUALITY_CONTROLLED_V2";
    row.observed_topology_mode = get_polar_topology_mode(config);
    row.expected_combined_hash = expected.combined_geometry_hash;
    row.actual_combined_hash = actual.combined_geometry_hash;
    row.exact_hash_match = row.expected_combined_hash == ...
        row.actual_combined_hash;
    row.bem_unknown_count = actual.total_panel_count;
    row.status = pass_fail(row.observed_topology_mode == ...
        "QUALITY_CONTROLLED_V2" && row.exact_hash_match);
    row.notes = "fresh L3 v5 mesh/operator preview; no matrix assembly";
    clear cleanup
end

function row = run_legacy_response(definition, suitePosition)
    omega = 1.5; % [rad/s]
    referenceTable = readtable(fullfile(definition.convergenceDirectory, ...
        'Phase3_1_Artifacts', 'Phase3_1_Formal_PFull_Frequency_Sweep.csv'), ...
        'TextType', 'string');
    selected = abs(referenceTable.omega_rad_s - omega) <= 64 * eps(omega);
    assert(nnz(selected) == 1, 'CRESTU:Phase32BLegacyResponseReference', ...
        'Frozen legacy 1.5 rad/s response row is missing.');
    reference = referenceTable(selected, :);
    runtimeDirectory = tempname(definition.artifactDirectory);
    mkdir(runtimeDirectory);
    cleanup = onCleanup(@() phase32.remove_runtime( ...
        definition, runtimeDirectory));
    [~, token] = fileparts(runtimeDirectory);
    caseName = sprintf('P32B_LegacyResponse_w1500_%s', token);
    bodyFile = fullfile(runtimeDirectory, [caseName, '_body.bmf']);
    configFile = fullfile(runtimeDirectory, [caseName, '.cfg']);
    baseline = read_config(definition.baselineConfig);
    generate_body_bmf(bodyFile, 2 * definition.sphereRadiusM, ...
        baseline.isx, baseline.isy, 7);
    write_legacy_response_config(configFile, caseName, bodyFile, omega, baseline);
    timer = tic;
    results = run_frequency_domain_case(configFile, ...
        'skipPhysicalDiagnostics', true);
    runtimeS = toc(timer);
    loaded = load(results.potential_cache_file, 'pot_cache');
    entry = loaded.pot_cache.entries(1);
    diagnostics = entry.solver_diagnostics;
    conditioning = cellfun(@(item) item.conditioning, diagnostics, ...
        'UniformOutput', false);
    rawRcond = min(cellfun(@(item) item.rawRcondEstimate, conditioning));
    scaledRcond = min(cellfun(@(item) item.scaledRcondEstimate, conditioning));
    maximumResidual = max([cellfun(@(item) ...
        item.radiationRelativeResidual, diagnostics); cellfun(@(item) ...
        item.diffractionRelativeResidual, diagnostics); cellfun(@(item) ...
        item.radiation.relativeResidualFrobeniusNorm, diagnostics); ...
        cellfun(@(item) item.diffraction.relativeResidualFrobeniusNorm, ...
        diagnostics)]);
    excitation = results.excitation(3, 1, 1);
    response = results.rao.complex(3, 1, 1);
    referenceExcitation = complex(reference.F3_real_N_per_m, ...
        reference.F3_imag_N_per_m);
    referenceResponse = complex(reference.RAO3_real_m_per_m, ...
        reference.RAO3_imag_m_per_m);
    differences = [relative_difference(results.added_mass(3, 3, 1), ...
        reference.A33_kg); ...
        relative_difference(results.damping(3, 3, 1), reference.B33_kg_s); ...
        relative_difference(excitation, referenceExcitation); ...
        relative_difference(response, referenceResponse)];
    geometryMatch = string(entry.mesh_audit.geometry.hash) == ...
        reference.combined_geometry_hash;
    cacheStatus = string(results.audit.frequencyEntries(1).cacheStatus);
    pass = geometryMatch && max(differences) <= 1.0e-10 && ...
        rawRcond >= definition.rawRcondMinimum && ...
        scaledRcond >= definition.scaledRcondMinimum && ...
        maximumResidual <= definition.relativeResidualMaximum && ...
        lower(cacheStatus) == "miss";
    row = empty_row();
    row.suite_position = suitePosition;
    row.test_id = "LEGACY_1P5_RESPONSE_REGRESSION";
    row.omega_rad_s = omega;
    row.topology_request = "OMITTED";
    row.observed_topology_mode = "LEGACY";
    row.expected_combined_hash = reference.combined_geometry_hash;
    row.actual_combined_hash = string(entry.mesh_audit.geometry.hash);
    row.exact_hash_match = geometryMatch;
    row.maximum_response_relative_difference = max(differences);
    row.raw_rcond = rawRcond;
    row.scaled_rcond = scaledRcond;
    row.maximum_relative_residual = maximumResidual;
    row.cache_status = cacheStatus;
    row.bem_unknown_count = entry.mesh_audit.bemUnknownCount;
    row.total_runtime_s = runtimeS;
    row.status = pass_fail(pass);
    row.notes = "fresh production solve compared with frozen Phase 3.1 row";
    clear cleanup
end

function write_q_only_config(filename, caseName, omega, qTop, qBottom, ...
        baseline, formal)
    write_config(filename, caseName, baseline.bodies(1).mesh_file, omega, ...
        qTop, qBottom, baseline, formal.wave.ndir, formal.wave.headings, false);
end

function write_legacy_response_config(filename, caseName, bodyFile, omega, baseline)
    write_config(filename, caseName, bodyFile, omega, 2.0, 1.5, ...
        baseline, 1, 0.0, true);
end

function write_config(filename, caseName, bodyFile, omega, qTop, qBottom, ...
        baseline, headingCount, headings, forceEnabled)
    fileID = fopen(filename, 'w');
    assert(fileID >= 0, 'CRESTU:Phase32BLegacyConfigOpen', ...
        'Cannot create legacy regression config.');
    cleanup = onCleanup(@() fclose(fileID));
    mass = baseline.mass_props(1);
    inertia = mass.inertia;
    inertiaTokens = [inertia(1, 1), inertia(2, 2), inertia(3, 3), ...
        -inertia(1, 2), -inertia(1, 3), -inertia(2, 3)];
    if forceEnabled
        forceFlags = '1 1 1 1 0';
    else
        forceFlags = '1 0 0 1 0';
    end
    fprintf(fileID, 'PARA1: Project Case Name\n%s\n\n', caseName);
    fprintf(fileID, 'PARA2: IPOTEN IFORCE IRAD IDIFF IDRIFT\n%s\n\n', forceFlags);
    fprintf(fileID, 'PARA3: FREQUENCY AND HEADING\n1 1 %.17g %.17g\n', omega, omega);
    fprintf(fileID, '%d\n', headingCount);
    fprintf(fileID, '%.17g ', headings);
    fprintf(fileID, '\n\nPARA4: ENVIRONMENT\n%.17g %.17g %.17g\n\n', ...
        baseline.water_depth, baseline.grav, baseline.rho);
    fprintf(fileID, 'PARA5: BODY\n1\n%d %s %.17g %.17g %.17g %.17g\n\n', ...
        baseline.bodies(1).id, bodyFile, baseline.bodies(1).pos, ...
        baseline.bodies(1).yaw);
    fprintf(fileID, 'PARA6: MASS PROPERTIES\n%.17g %.17g %.17g %.17g\n', ...
        mass.mass, mass.cg);
    fprintf(fileID, '%.17g %.17g %.17g %.17g %.17g %.17g\n\n', inertiaTokens);
    fprintf(fileID, 'PARA7: MODES\n');
    fprintf(fileID, '%d ', baseline.calc_modes);
    fprintf(fileID, '\n\nPARA8: FULL DOMAIN\n%d %d\n\n', baseline.isx, baseline.isy);
    fprintf(fileID, ['PARA9: DOMAIN\n%.17g %d %d %.17g %.17g ', ...
        '%.17g %d\n\n'], baseline.z_tol, baseline.fs.nr_near, ...
        baseline.fs.nr_sponge, baseline.fs.r_inner, baseline.fs.r_outer, ...
        baseline.fs.sponge_ratio, baseline.fs.nz_farfield);
    fprintf(fileID, 'PARA10: SPONGE\n%.17g\n\n', baseline.fs.mu0);
    fprintf(fileID, ['PARA11: FIXED OUTER-DOMAIN WATERLINE\n', ...
        'OUTER_WATERLINE_MESH %s\n\n'], bodyFile);
    fprintf(fileID, ['PARA12: PHASE2.2 OUTER TRUNCATION AUDIT\n', ...
        'OUTER_TRUNCATION_Q %.17g %.17g\n\n'], qTop, qBottom);
    fprintf(fileID, ['PARA13: LEGACY REGRESSION DIAGNOSTICS\n', ...
        'RAW_DIAGNOSTICS 1\nPERSIST_FULL_PHI 0\n']);
    clear cleanup
end

function value = relative_difference(actual, expected)
    value = abs(actual - expected) / max([abs(actual), abs(expected), eps]);
end

function row = empty_row()
    row = struct('suite_position', "", 'test_id', "", ...
        'omega_rad_s', NaN, 'topology_request', "", ...
        'observed_topology_mode', "", 'expected_combined_hash', "", ...
        'actual_combined_hash', "", 'exact_hash_match', false, ...
        'maximum_response_relative_difference', NaN, ...
        'raw_rcond', NaN, 'scaled_rcond', NaN, ...
        'maximum_relative_residual', NaN, 'cache_status', "NOT_RUN", ...
        'bem_unknown_count', 0, 'total_runtime_s', NaN, ...
        'status', "NOT_RUN", 'notes', "");
end

function text = pass_fail(value)
    if value
        text = "PASS";
    else
        text = "FAIL";
    end
end
