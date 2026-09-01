function summary = Finalize_Phase3_2B_Resource_Plan( ...
        peakWorkingSetBytes, processWallClockS)
% FINALIZE_PHASE3_2B_RESOURCE_PLAN Bind measured process resources to plan.

    validateattributes(peakWorkingSetBytes, {'numeric'}, ...
        {'scalar', 'real', 'nonnegative', 'finite'});
    validateattributes(processWallClockS, {'numeric'}, ...
        {'scalar', 'real', 'nonnegative', 'finite'});
    definition = phase32b.initialize();
    planFile = fullfile(definition.artifactDirectory, ...
        'Phase3_2B_Resource_Plan.csv');
    pilotFile = fullfile(definition.artifactDirectory, ...
        'Phase3_2B_Resource_Pilot.csv');
    plan = readtable(planFile, 'TextType', 'string');
    pilot = readtable(pilotFile, 'TextType', 'string');
    assert(height(pilot) == 1 && nnz(plan.maximum_N_case) == 1, ...
        'CRESTU:Phase32BResourcePilotShape', ...
        'One maximum-N pilot record is required.');
    selected = logical(plan.maximum_N_case);
    plan.measured_peak_working_set_bytes(selected) = peakWorkingSetBytes;
    plan.measured_process_wall_clock_s(selected) = processWallClockS;
    plan.pilot_solve_status(selected) = pilot.solve_status;
    if pilot.solve_status == "PASS"
        gate = "PASS_SERIAL_L3";
    elseif pilot.solve_status == "BLOCKED_RESOURCE"
        gate = "BLOCKED_RESOURCE";
    else
        gate = "FAIL_NONRESOURCE";
    end
    plan.resource_gate_status(:) = gate;
    memoryFraction = peakWorkingSetBytes / ...
        plan.total_physical_memory_bytes(find(selected, 1));
    if isfinite(memoryFraction) && memoryFraction > 0.70
        policy = "FORMAL_DEGRADE_2_FINEST_MEMORY_GT_70_PERCENT";
    else
        policy = "FULL_FORMAL_PENDING_SENTINEL_RUNTIME_FORECAST";
    end
    plan.formal_resource_policy(:) = policy;
    writetable(plan, planFile);
    summary = struct('schemaVersion', 1, 'resourceGateStatus', char(gate), ...
        'peakWorkingSetBytes', peakWorkingSetBytes, ...
        'processWallClockS', processWallClockS, ...
        'memoryFraction', memoryFraction, ...
        'formalResourcePolicy', char(policy));
    fprintf('[RESULT] RESOURCE_GATE = %s\n', gate);
    fprintf('[RESULT] PILOT_PEAK_WORKING_SET_GIB = %.6g\n', ...
        peakWorkingSetBytes / 2^30);
end
