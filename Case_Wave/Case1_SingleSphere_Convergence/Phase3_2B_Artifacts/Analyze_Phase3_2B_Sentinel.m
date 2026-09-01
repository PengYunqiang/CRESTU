function summary = Analyze_Phase3_2B_Sentinel()
% ANALYZE_PHASE3_2B_SENTINEL Apply frozen adjacent-mesh response gates.

    %% Stage 1: Validate and enrich the 18 raw clean-solve records

    definition = phase32b.initialize();
    responseFile = fullfile(definition.artifactDirectory, ...
        'Phase3_2B_Sentinel_Response.csv');
    response = readtable(responseFile, 'TextType', 'string');
    assert(height(response) == 18 && ...
        numel(unique(response.case_id)) == 18, ...
        'CRESTU:Phase32BSentinelCaseCount', ...
        'Sentinel evidence must contain 18 unique case IDs.');
    response = sortrows(response, {'mesh_level_index', 'omega_rad_s'});
    response = add_unwrapped_phases(response);
    writetable(response, responseFile);

    %% Stage 2: Export conditioning/residual and measured mesh-cost tables

    maximumResidual = max([response.radiation_linear_relative_residual, ...
        response.diffraction_linear_relative_residual, ...
        response.radiation_algebraic_relative_residual, ...
        response.diffraction_algebraic_relative_residual], [], 2);
    condition = response(:, {'case_id', 'mesh_level', 'mesh_level_index', ...
        'omega_rad_s', 'bem_unknown_count', 'matrix_rows', 'matrix_columns', ...
        'raw_rcond', 'scaled_rcond', ...
        'radiation_linear_relative_residual', ...
        'diffraction_linear_relative_residual', ...
        'radiation_algebraic_relative_residual', ...
        'diffraction_algebraic_relative_residual', ...
        'maximum_abs_outer_ntheta', ...
        'maximum_outer_column_residual', 'all_finite', ...
        'solver_warning_count', 'cache_status', 'algebraic_gate_status'});
    condition.maximum_relative_residual = maximumResidual;
    writetable(condition, fullfile(definition.artifactDirectory, ...
        'Phase3_2B_Condition_Residual.csv'));
    meshCost = build_mesh_cost(definition, response);
    writetable(meshCost, fullfile(definition.artifactDirectory, ...
        'Phase3_2B_Mesh_Cost.csv'));

    %% Stage 3: Compute scalar, component, and complex adjacent differences

    frequencies = sort(unique(response.omega_rad_s));
    quantityDefinitions = build_quantity_definitions(definition);
    rows = repmat(empty_convergence_row(), 0, 1);
    for frequencyIndex = 1:numel(frequencies)
        omega = frequencies(frequencyIndex); % [rad/s]
        frequencyRows = response(abs(response.omega_rad_s - omega) <= ...
            64 * eps(max(1.0, omega)), :);
        frequencyRows = sortrows(frequencyRows, 'mesh_level_index');
        assert(height(frequencyRows) == 3 && ...
            isequal(frequencyRows.mesh_level_index, (1:3)'), ...
            'CRESTU:Phase32BSentinelLevels', ...
            'Each sentinel frequency requires L1/L2/L3 exactly once.');
        for quantityIndex = 1:numel(quantityDefinitions)
            definitionRow = quantityDefinitions(quantityIndex);
            row = analyze_quantity( ...
                frequencyRows, definitionRow, definition, omega);
            rows(end + 1, 1) = row; %#ok<AGROW>
        end
    end
    convergence = struct2table(rows);
    writetable(convergence, fullfile(definition.artifactDirectory, ...
        'Phase3_2B_Sentinel_Convergence.csv'));

    %% Stage 4: Combine frozen algebraic and response continuation gates

    algebraicPass = all(response.algebraic_gate_status == "PASS") && ...
        all(response.preassembly_status == "PASS") && ...
        all(response.v5_mesh_reproducibility == "PASS") && ...
        all(lower(response.cache_status) == "miss");
    mandatory = convergence.gate_role == "MANDATORY";
    responseFail = any(convergence.gate_status(mandatory) == "FAIL");
    responseBlocked = any(contains(response.solve_status, "BLOCKED"));
    resonanceDeferred = any(convergence.gate_status(mandatory) == ...
        "DEFER_RESONANCE_TRACKING");
    if responseBlocked
        responseGate = "BLOCKED";
    elseif responseFail
        responseGate = "FAIL";
    elseif resonanceDeferred
        responseGate = "PASS_WITH_RESONANCE_DEFERRED";
    else
        responseGate = "PASS";
    end
    update_resource_policy(definition, responseGate);
    write_figures(definition, response, convergence);
    write_report(definition, response, convergence, algebraicPass, responseGate);
    summary = struct('schemaVersion', 1, ...
        'caseCount', height(response), ...
        'cleanSolveCount', nnz(lower(response.cache_status) == "miss"), ...
        'algebraicPassCount', ...
            nnz(response.algebraic_gate_status == "PASS"), ...
        'sentinelAlgebraicGate', pass_fail(algebraicPass), ...
        'sentinelResponseGate', char(responseGate), ...
        'mandatoryComparisonCount', nnz(mandatory), ...
        'mandatoryFailCount', ...
            nnz(convergence.gate_status(mandatory) == "FAIL"), ...
        'resonanceDeferredCount', ...
            nnz(convergence.gate_status(mandatory) == ...
            "DEFER_RESONANCE_TRACKING"));
    fprintf('[RESULT] SENTINEL_ALGEBRAIC_GATE = %s\n', ...
        summary.sentinelAlgebraicGate);
    fprintf('[RESULT] SENTINEL_RESPONSE_GATE = %s\n', responseGate);
end

function response = add_unwrapped_phases(response)
    response.F3_phase_unwrapped_deg = NaN(height(response), 1);
    response.RAO3_phase_unwrapped_deg = NaN(height(response), 1);
    response.D3_phase_unwrapped_deg = NaN(height(response), 1);
    for levelIndex = 1:3
        selected = response.mesh_level_index == levelIndex;
        [~, order] = sort(response.omega_rad_s(selected));
        indices = find(selected);
        indices = indices(order);
        response.F3_phase_unwrapped_deg(indices) = rad2deg(unwrap(deg2rad( ...
            response.F3_phase_deg(indices))));
        response.RAO3_phase_unwrapped_deg(indices) = rad2deg(unwrap(deg2rad( ...
            response.RAO3_phase_deg(indices))));
        response.D3_phase_unwrapped_deg(indices) = rad2deg(unwrap(deg2rad( ...
            response.D3_phase_deg(indices))));
    end
end

function definitions = build_quantity_definitions(definition)
    template = struct('name', "", 'kind', "", 'gateRole', "", ...
        'floor', NaN, 'gate', NaN, 'realField', "", 'imagField', "", ...
        'phaseField', "");
    definitions = repmat(template, 11, 1);
    definitions(1) = scalar_definition("A33", "MANDATORY", ...
        definition.responseFloors.A33Kg, ...
        definition.responseRelativeGates.A33, "A33_kg");
    definitions(2) = scalar_definition("B33", "MANDATORY", ...
        definition.responseFloors.B33KgPerS, ...
        definition.responseRelativeGates.B33, "B33_kg_s");
    definitions(3) = complex_definition("F3", "MANDATORY", ...
        definition.responseFloors.F3NPerM, ...
        definition.responseRelativeGates.F3, ...
        "F3_real_N_per_m", "F3_imag_N_per_m", ...
        "F3_phase_unwrapped_deg");
    definitions(4) = component_definition("F3_REAL", ...
        definition.responseFloors.F3NPerM, "F3_real_N_per_m");
    definitions(5) = component_definition("F3_IMAG", ...
        definition.responseFloors.F3NPerM, "F3_imag_N_per_m");
    definitions(6) = complex_definition("RAO3", "MANDATORY", ...
        definition.responseFloors.RAO3MPerM, ...
        definition.responseRelativeGates.RAO3, ...
        "RAO3_real_m_per_m", "RAO3_imag_m_per_m", ...
        "RAO3_phase_unwrapped_deg");
    definitions(7) = component_definition("RAO3_REAL", ...
        definition.responseFloors.RAO3MPerM, "RAO3_real_m_per_m");
    definitions(8) = component_definition("RAO3_IMAG", ...
        definition.responseFloors.RAO3MPerM, "RAO3_imag_m_per_m");
    definitions(9) = complex_definition("D3", "DIAGNOSTIC", ...
        definition.responseFloors.D3NPerM, NaN, ...
        "D3_real_N_per_m", "D3_imag_N_per_m", ...
        "D3_phase_unwrapped_deg");
    definitions(10) = component_definition("D3_REAL", ...
        definition.responseFloors.D3NPerM, "D3_real_N_per_m");
    definitions(11) = component_definition("D3_IMAG", ...
        definition.responseFloors.D3NPerM, "D3_imag_N_per_m");
end

function value = scalar_definition(name, role, floorValue, gate, field)
    value = struct('name', name, 'kind', "REAL_SCALAR", ...
        'gateRole', role, 'floor', floorValue, 'gate', gate, ...
        'realField', field, 'imagField', "", 'phaseField', "");
end

function value = component_definition(name, floorValue, field)
    value = scalar_definition(name, "DIAGNOSTIC_GCI_ONLY", ...
        floorValue, NaN, field);
end

function value = complex_definition(name, role, floorValue, gate, ...
        realField, imagField, phaseField)
    value = struct('name', name, 'kind', "COMPLEX", 'gateRole', role, ...
        'floor', floorValue, 'gate', gate, 'realField', realField, ...
        'imagField', imagField, 'phaseField', phaseField);
end

function row = analyze_quantity(response, quantity, definition, omega)
    if quantity.kind == "COMPLEX"
        values = complex(response.(quantity.realField), ...
            response.(quantity.imagField));
        phases = response.(quantity.phaseField);
    else
        values = response.(quantity.realField);
        phases = NaN(3, 1);
    end
    delta12 = values(2) - values(1);
    delta23 = values(3) - values(2);
    absolute12 = abs(delta12);
    absolute23 = abs(delta23);
    relative12 = absolute12 / max(abs(values(2)), quantity.floor);
    relative23 = absolute23 / max(abs(values(3)), quantity.floor);
    relative23LowFloor = absolute23 / ...
        max(abs(values(3)), 0.1 * quantity.floor);
    relative23HighFloor = absolute23 / ...
        max(abs(values(3)), 10.0 * quantity.floor);
    if absolute12 <= quantity.floor * 1.0e-6
        convergenceRatio = NaN;
    else
        convergenceRatio = absolute23 / absolute12;
    end

    [monotonicity, trajectory, turnAngle] = ...
        trajectory_status(values, quantity.kind, delta12, delta23);
    [observedOrder, richardson, gci, asymptotic] = ...
        scalar_extrapolation(values, quantity, definition);
    [classification, gateStatus] = classify_result(quantity, omega, ...
        relative23, absolute12, absolute23, monotonicity, trajectory, ...
        observedOrder, gci);
    if isfinite(quantity.gate)
        decisions = [relative23 <= quantity.gate, ...
            relative23LowFloor <= quantity.gate, ...
            relative23HighFloor <= quantity.gate];
        floorSensitivity = stable_status(all(decisions == decisions(1)));
    else
        floorSensitivity = "DIAGNOSTIC_NO_GATE";
    end

    row = empty_convergence_row();
    row.omega_rad_s = omega;
    row.quantity = quantity.name;
    row.value_kind = quantity.kind;
    row.gate_role = quantity.gateRole;
    row.Q_floor = quantity.floor;
    row.relative_gate = quantity.gate;
    row.L1_real = real(values(1)); row.L1_imag = imag(values(1));
    row.L1_magnitude = abs(values(1)); row.L1_phase_unwrapped_deg = phases(1);
    row.L2_real = real(values(2)); row.L2_imag = imag(values(2));
    row.L2_magnitude = abs(values(2)); row.L2_phase_unwrapped_deg = phases(2);
    row.L3_real = real(values(3)); row.L3_imag = imag(values(3));
    row.L3_magnitude = abs(values(3)); row.L3_phase_unwrapped_deg = phases(3);
    row.L1_L2_absolute_difference = absolute12;
    row.L1_L2_relative_difference = relative12;
    row.L2_L3_absolute_difference = absolute23;
    row.L2_L3_relative_difference = relative23;
    row.L2_L3_relative_difference_floor_x0p1 = relative23LowFloor;
    row.L2_L3_relative_difference_floor_x10 = relative23HighFloor;
    row.floor_sensitivity_status = floorSensitivity;
    row.convergence_ratio = convergenceRatio;
    row.monotonicity = monotonicity;
    row.complex_trajectory = trajectory;
    row.complex_turn_angle_deg = turnAngle;
    row.observed_order = observedOrder;
    row.richardson_extrapolated = richardson;
    row.fine_grid_GCI = gci;
    row.asymptotic_range = asymptotic;
    row.classification = classification;
    row.gate_status = gateStatus;
end

function [monotonicity, trajectory, turnAngle] = ...
        trajectory_status(values, kind, delta12, delta23)
    if kind == "COMPLEX"
        if abs(delta12) == 0 || abs(delta23) == 0
            trajectory = "DEGENERATE";
            turnAngle = NaN;
        else
            turnAngle = abs(rad2deg(angle(delta23 / delta12)));
            if turnAngle <= 90
                trajectory = "DIRECT";
            else
                trajectory = "REVERSING";
            end
        end
        monotonicity = "NOT_APPLICABLE_COMPLEX";
    else
        trajectory = "NOT_APPLICABLE_REAL";
        turnAngle = NaN;
        products = (values(2) - values(1)) * (values(3) - values(2));
        if products > 0
            monotonicity = "MONOTONIC";
        elseif products < 0
            monotonicity = "OSCILLATORY";
        else
            monotonicity = "PLATEAU_OR_DEGENERATE";
        end
    end
end

function [order, extrapolated, gci, status] = ...
        scalar_extrapolation(values, quantity, definition)
    order = NaN; extrapolated = NaN; gci = NaN;
    status = "NOT_DEMONSTRATED";
    if quantity.kind == "COMPLEX" || any(~isfinite(values)) || ...
            any(values == 0) || min(values) < 0 && max(values) > 0
        return
    end
    delta12 = values(1) - values(2);
    delta23 = values(2) - values(3);
    if delta12 * delta23 <= 0 || delta23 == 0
        return
    end
    h = [definition.levels.eta];
    r21 = h(1) / h(2);
    r32 = h(2) / h(3);
    target = abs(delta12 / delta23);
    equation = @(p) r32 .^ p .* (r21 .^ p - 1) ./ ...
        (r32 .^ p - 1) - target;
    lowerValue = equation(0.1);
    upperValue = equation(10.0);
    if ~isfinite(lowerValue) || ~isfinite(upperValue) || ...
            lowerValue * upperValue > 0
        return
    end
    order = fzero(equation, [0.1, 10.0]);
    if ~isfinite(order) || order < 0.1 || order > 10.0
        order = NaN;
        return
    end
    extrapolated = values(3) + (values(3) - values(2)) / ...
        (r32 ^ order - 1);
    gci = 1.25 * abs(values(3) - values(2)) / ...
        max(abs(values(3)), quantity.floor) / (r32 ^ order - 1);
    status = "DEMONSTRATED_FOR_REAL_SCALAR";
end

function [classification, gateStatus] = classify_result(quantity, omega, ...
        relative23, absolute12, absolute23, monotonicity, trajectory, ...
        observedOrder, gci)
    if ~isfinite(quantity.gate)
        classification = "DIAGNOSTIC_ONLY";
        gateStatus = "NOT_APPLICABLE";
        return
    end
    if relative23 <= quantity.gate
        gateStatus = "PASS";
        if quantity.kind == "COMPLEX" && trajectory == "REVERSING" || ...
                quantity.kind ~= "COMPLEX" && monotonicity == "OSCILLATORY"
            classification = "OSCILLATORY_CONVERGENCE";
        elseif isfinite(observedOrder) && isfinite(gci) && gci <= quantity.gate
            classification = "CONVERGED";
        else
            classification = "CONVERGED_WITHOUT_ASYMPTOTIC_PROOF";
        end
        return
    end
    resonanceCandidate = quantity.name == "RAO3" && ...
        any(abs(omega - [1.4, 1.5]) <= 64 * eps(max(1.0, omega))) && ...
        relative23 <= 0.25 && absolute23 < absolute12 && ...
        trajectory == "DIRECT";
    if resonanceCandidate
        classification = "RESONANCE_SHIFT";
        gateStatus = "DEFER_RESONANCE_TRACKING";
    elseif absolute23 < absolute12
        classification = "PRE_ASYMPTOTIC";
        gateStatus = "FAIL";
    else
        classification = "NONCONVERGENT";
        gateStatus = "FAIL";
    end
end

function status = stable_status(value)
    if value
        status = "STABLE_X0P1_TO_X10";
    else
        status = "FLOOR_SENSITIVE";
    end
end

function meshCost = build_mesh_cost(definition, response)
    resource = readtable(fullfile(definition.artifactDirectory, ...
        'Phase3_2B_Resource_Plan.csv'), 'TextType', 'string');
    pilot = resource(logical(resource.maximum_N_case), :);
    meshCost = response(:, {'case_id', 'mesh_level', 'mesh_level_index', ...
        'omega_rad_s', 'bem_unknown_count', 'assembly_runtime_s', ...
        'total_runtime_s'});
    meshCost.one_dense_matrix_bytes = ...
        16 .* meshCost.bem_unknown_count .^ 2;
    meshCost.A_L_U_lower_bound_bytes = ...
        3 .* meshCost.one_dense_matrix_bytes;
    meshCost.peak_memory_estimated_from_pilot_bytes = ...
        pilot.measured_peak_working_set_bytes .* ...
        (meshCost.bem_unknown_count ./ pilot.bem_unknown_count) .^ 2;
    meshCost.peak_memory_basis = repmat( ...
        "N2_CALIBRATED_FROM_MAXIMUM_N_PROCESS_PEAK", height(meshCost), 1);
    meshCost.resource_pilot_peak_working_set_bytes = repmat( ...
        pilot.measured_peak_working_set_bytes, height(meshCost), 1);
    meshCost.maximum_concurrent_L3_BEM_solves = repmat( ...
        definition.maximumConcurrentL3Solves, height(meshCost), 1);
end

function write_figures(definition, response, convergence)
    levels = unique(response.mesh_level, 'stable');
    colors = lines(numel(levels));
    figureHandle = figure('Visible', 'off', 'Color', 'w', ...
        'Position', [100, 100, 1200, 850]);
    quantities = {'A33_kg', 'B33_kg_s', 'F3_abs_N_per_m', ...
        'RAO3_abs_m_per_m'};
    labels = {'A33 [kg]', 'B33 [kg/s]', '|F3| [N/m]', '|RAO3| [m/m]'};
    for plotIndex = 1:4
        subplot(2, 2, plotIndex); hold on; grid on;
        for levelIndex = 1:numel(levels)
            selected = response.mesh_level == levels(levelIndex);
            data = sortrows(response(selected, :), 'omega_rad_s');
            plot(data.omega_rad_s, data.(quantities{plotIndex}), '-o', ...
                'LineWidth', 1.4, 'Color', colors(levelIndex, :), ...
                'DisplayName', levels(levelIndex));
        end
        xlabel('\omega [rad/s]'); ylabel(labels{plotIndex}); legend('Location', 'best');
    end
    exportgraphics(figureHandle, fullfile(definition.artifactDirectory, ...
        'Phase3_2B_Sentinel_Response_Curves.png'), 'Resolution', 180);
    close(figureHandle);

    figureHandle = figure('Visible', 'off', 'Color', 'w', ...
        'Position', [100, 100, 1000, 650]); hold on; grid on;
    mandatory = convergence.gate_role == "MANDATORY";
    names = unique(convergence.quantity(mandatory), 'stable');
    for nameIndex = 1:numel(names)
        selected = mandatory & convergence.quantity == names(nameIndex);
        data = sortrows(convergence(selected, :), 'omega_rad_s');
        plot(data.omega_rad_s, 100 * data.L2_L3_relative_difference, ...
            '-o', 'LineWidth', 1.4, 'DisplayName', names(nameIndex));
    end
    xlabel('\omega [rad/s]'); ylabel('L2 to L3 relative difference [%]');
    legend('Location', 'best');
    exportgraphics(figureHandle, fullfile(definition.artifactDirectory, ...
        'Phase3_2B_Sentinel_Adjacent_Sensitivity.png'), 'Resolution', 180);
    close(figureHandle);
end

function write_report(definition, response, convergence, algebraicPass, responseGate)
    fileID = fopen(fullfile(definition.artifactDirectory, ...
        'PHASE3_2B_RESPONSE_CONVERGENCE_CN.md'), 'w');
    assert(fileID >= 0, 'CRESTU:Phase32BResponseReportOpen', ...
        'Cannot create the response convergence report.');
    cleanup = onCleanup(@() fclose(fileID));
    mandatory = convergence.gate_role == "MANDATORY";
    fprintf(fileID, '# CRESTU-1F Phase 3.2B 响应收敛\n\n');
    fprintf(fileID, '```ini\nV5_MESH_REPRODUCIBILITY = %s\n', ...
        pass_fail(all(response.v5_mesh_reproducibility == "PASS")));
    fprintf(fileID, 'SENTINEL_CLEAN_SOLVES = %d/18\n', ...
        nnz(lower(response.cache_status) == "miss"));
    fprintf(fileID, 'SENTINEL_ALGEBRAIC_GATE = %s\n', pass_fail(algebraicPass));
    fprintf(fileID, 'SENTINEL_RESPONSE_GATE = %s\n', responseGate);
    fprintf(fileID, 'MANDATORY_COMPARISONS = %d\n', nnz(mandatory));
    fprintf(fileID, 'MANDATORY_FAILS = %d\n', ...
        nnz(convergence.gate_status(mandatory) == "FAIL"));
    if any(responseGate == ["FAIL", "BLOCKED"])
        fprintf(fileID, ['RESONANCE_CONVERGENCE = ', ...
            'NOT_RUN_BLOCKED_SENTINEL_RESPONSE_%s\n'], responseGate);
        fprintf(fileID, ['FORMAL_MESH_CONVERGENCE = ', ...
            'NOT_RUN_BLOCKED_SENTINEL_RESPONSE_%s\n'], responseGate);
    else
        fprintf(fileID, 'RESONANCE_CONVERGENCE = PENDING\n');
        fprintf(fileID, 'FORMAL_MESH_CONVERGENCE = PENDING\n');
    end
    fprintf(fileID, 'PHYSICAL_VALIDATION = PARTIAL\n');
    fprintf(fileID, 'READY_FOR_PRODUCTION = NO\n```\n\n');
    fprintf(fileID, ['全部 sentinel 均为唯一临时目录中的真实 clean solve；', ...
        'BMF/MAT/cache 在提取 compact CSV 后删除。复数 F3/RAO3 以 complex norm ', ...
        '判定，Re/Im 的 GCI 仅作为独立实分量诊断。\n\n']);
    fprintf(fileID, '分类计数：\n\n');
    names = unique(convergence.classification);
    for index = 1:numel(names)
        fprintf(fileID, '- `%s`: %d\n', names(index), ...
            nnz(convergence.classification == names(index)));
    end
    clear cleanup
end

function update_resource_policy(definition, responseGate)
    resourceFile = fullfile(definition.artifactDirectory, ...
        'Phase3_2B_Resource_Plan.csv');
    resource = readtable(resourceFile, 'TextType', 'string');
    if responseGate == "FAIL"
        resource.formal_resource_policy(:) = ...
            "NOT_RUN_BLOCKED_SENTINEL_RESPONSE_FAIL";
    elseif responseGate == "BLOCKED"
        resource.formal_resource_policy(:) = ...
            "NOT_RUN_BLOCKED_SENTINEL_RESPONSE_BLOCKED";
    else
        resource.formal_resource_policy(:) = ...
            "FULL_FORMAL_PENDING_RESONANCE_AND_REVIEWER_GATES";
    end
    writetable(resource, resourceFile);
end

function row = empty_convergence_row()
    row = struct('omega_rad_s', NaN, 'quantity', "", ...
        'value_kind', "", 'gate_role', "", 'Q_floor', NaN, ...
        'relative_gate', NaN, 'L1_real', NaN, 'L1_imag', NaN, ...
        'L1_magnitude', NaN, 'L1_phase_unwrapped_deg', NaN, ...
        'L2_real', NaN, 'L2_imag', NaN, 'L2_magnitude', NaN, ...
        'L2_phase_unwrapped_deg', NaN, 'L3_real', NaN, ...
        'L3_imag', NaN, 'L3_magnitude', NaN, ...
        'L3_phase_unwrapped_deg', NaN, ...
        'L1_L2_absolute_difference', NaN, ...
        'L1_L2_relative_difference', NaN, ...
        'L2_L3_absolute_difference', NaN, ...
        'L2_L3_relative_difference', NaN, ...
        'L2_L3_relative_difference_floor_x0p1', NaN, ...
        'L2_L3_relative_difference_floor_x10', NaN, ...
        'floor_sensitivity_status', "", 'convergence_ratio', NaN, ...
        'monotonicity', "", 'complex_trajectory', "", ...
        'complex_turn_angle_deg', NaN, 'observed_order', NaN, ...
        'richardson_extrapolated', NaN, 'fine_grid_GCI', NaN, ...
        'asymptotic_range', "NOT_DEMONSTRATED", ...
        'classification', "", 'gate_status', "");
end

function text = pass_fail(value)
    if value
        text = "PASS";
    else
        text = "FAIL";
    end
end
