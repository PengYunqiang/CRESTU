function plan = Prepare_Phase3_2B_Resource_Plan()
% PREPARE_PHASE3_2B_RESOURCE_PLAN Estimate dense solve cost before batching.

    %% Stage 1: Read actual frozen N and the accepted bridge timing

    definition = phase32b.initialize();
    frozen = readtable(fullfile(definition.artifactDirectory, ...
        'Phase3_2B_V5_Frozen_Mesh_Manifest.csv'), 'TextType', 'string');
    bridge = readtable(fullfile(definition.phase32aDirectory, ...
        'Phase3_2A_BEM_Bridge_Smoke.csv'), 'TextType', 'string');
    bridgeRow = bridge(bridge.topology_mode == "QUALITY_CONTROLLED_V2", :);
    assert(height(bridgeRow) == 1 && bridgeRow.solve_status == "PASS", ...
        'CRESTU:Phase32BResourceBridge', ...
        'The accepted V2 bridge resource reference is missing.');
    totalPhysicalMemoryBytes = query_total_physical_memory();

    %% Stage 2: Calculate transparent matrix-memory lower bounds

    plan = table();
    plan.mesh_level = frozen.mesh_level;
    plan.mesh_level_index = frozen.mesh_level_index;
    plan.omega_rad_s = frozen.omega_rad_s;
    plan.bem_unknown_count = frozen.total_panel_count;
    plan.matrix_type = repmat("DENSE_COMPLEX_DOUBLE", height(frozen), 1);
    plan.bytes_per_complex_entry = repmat(16, height(frozen), 1);
    plan.one_dense_matrix_bytes = 16 .* frozen.total_panel_count .^ 2;
    plan.one_dense_matrix_gib = plan.one_dense_matrix_bytes ./ 2^30;
    plan.A_L_U_lower_bound_bytes = 3 .* plan.one_dense_matrix_bytes;
    plan.A_L_U_lower_bound_gib = plan.A_L_U_lower_bound_bytes ./ 2^30;
    plan.conditioning_five_matrix_bound_bytes = ...
        5 .* plan.one_dense_matrix_bytes;
    plan.conditioning_five_matrix_bound_gib = ...
        plan.conditioning_five_matrix_bound_bytes ./ 2^30;
    plan.bridge_unknown_count = repmat( ...
        bridgeRow.bem_unknown_count, height(frozen), 1);
    plan.bridge_assembly_runtime_s = repmat( ...
        bridgeRow.assembly_runtime_s, height(frozen), 1);
    plan.bridge_total_runtime_s = repmat( ...
        bridgeRow.total_runtime_s, height(frozen), 1);
    plan.bridge_peak_memory_status = repmat( ...
        "NOT_RECORDED_PILOT_REQUIRED", height(frozen), 1);
    plan.estimated_assembly_runtime_s_N2 = bridgeRow.assembly_runtime_s .* ...
        (frozen.total_panel_count ./ bridgeRow.bem_unknown_count) .^ 2;
    plan.estimated_total_runtime_s_N3 = bridgeRow.total_runtime_s .* ...
        (frozen.total_panel_count ./ bridgeRow.bem_unknown_count) .^ 3;
    plan.total_physical_memory_bytes = repmat( ...
        totalPhysicalMemoryBytes, height(frozen), 1);
    maximumN = max(frozen.total_panel_count);
    plan.maximum_N_case = frozen.total_panel_count == maximumN;
    plan.resource_pilot_required = plan.maximum_N_case;
    plan.maximum_concurrent_L3_BEM_solves = repmat( ...
        definition.maximumConcurrentL3Solves, height(frozen), 1);
    plan.measured_peak_working_set_bytes = NaN(height(frozen), 1);
    plan.measured_process_wall_clock_s = NaN(height(frozen), 1);
    plan.pilot_solve_status = repmat("PENDING", height(frozen), 1);
    plan.resource_gate_status = repmat("PENDING_PILOT", height(frozen), 1);
    plan.formal_resource_policy = repmat( ...
        "PENDING_SENTINEL_MEASUREMENTS", height(frozen), 1);
    writetable(plan, fullfile(definition.artifactDirectory, ...
        'Phase3_2B_Resource_Plan.csv'));
    fprintf('[RESULT] MAXIMUM_N = %d at omega=%.6g, %s\n', maximumN, ...
        frozen.omega_rad_s(plan.maximum_N_case), ...
        frozen.mesh_level(plan.maximum_N_case));
    fprintf('[RESULT] MAX_CONCURRENT_L3_BEM_SOLVES = %d\n', ...
        definition.maximumConcurrentL3Solves);
end

function bytes = query_total_physical_memory()
    bytes = NaN;
    try
        [~, systemView] = memory;
        bytes = systemView.PhysicalMemory.Total;
    catch
        fprintf('[WARN] MATLAB physical-memory query unavailable.\n');
    end
end
