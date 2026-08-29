function afterTable = Run_Phase2_2_Orientation_After()
% RUN_PHASE2_2_ORIENTATION_AFTER Verify the atomic double-layer repair.
%
% The mixed cube uses the production orientation helper. Only T2 is
% expected to pass after this change; all other analytic defects must remain
% visible and failing until separately authorized.

    testDirectory = fileparts(mfilename('fullpath'));
    artifactDirectory = fileparts(testDirectory);
    convergenceDirectory = fileparts(artifactDirectory);
    caseWaveDirectory = fileparts(convergenceDirectory);
    repositoryRoot = fileparts(caseWaveDirectory);
    evidenceDirectory = fullfile(testDirectory, 'Evidence');
    addpath(fullfile(repositoryRoot, 'Source Code', '2.Mesh'));
    addpath(fullfile(repositoryRoot, 'Source Code', '4.Potential'));
    addpath(fullfile(repositoryRoot, 'Source Code', '5.Force'));

    beforeTable = Run_Phase2_2_Analytic_FailBefore();
    panels = cube_panels_current_mixed(1.0);
    [normals, centers] = panel_normals_and_centers(panels);
    stats = struct('total_body_panels', 0, 'fs_panels', 1, ...
        'seabed_panels', 1, 'farfield_panels', 4);
    orientation = get_rankine_source_orientation(stats, normals, panels);
    canonicalD = assemble_canonical_double_layer( ...
        panels, normals, centers, orientation.columnSigns);
    repairedResidual = norm(canonicalD * ones(6, 1), inf);

    unaffectedIds = ["T3_FS_FINITE_DEPTH_NU", ...
        "T5_CYLINDER_STORED_DECAY", ...
        "T5_CYLINDER_GREEN_OUTWARD_SIGN", ...
        "T5_FRUSTUM_STORED_PROJECTION", ...
        "T5_COMPLEX_KAPPA_INTERFACE", ...
        "T6_GREEN_RHS_ORIENTATION", ...
        "T10_PRESSURE_TRACTION_SIGN", ...
        "T7_TRANSLATING_SPHERE_ADDED_MASS"];
    rows = repmat(empty_after_result(), 3 + numel(unaffectedIds), 1);
    rows(1) = make_after_result('T2_MIXED_CANONICAL_PRODUCTION_HELPER', ...
        repairedResidual, 1e-10, true, 'REPAIRED_PASS', ...
        'Production helper maps only FS dGdn by -1; D*1 returns to machine zero.');
    signError = norm(orientation.columnSigns - [-1; 1; 1; 1; 1; 1], inf);
    rows(2) = make_after_result('T2_COMPONENT_SIGN_CONTRACT', ...
        signError, 0.0, true, 'REPAIRED_PASS', ...
        'Component contract is body +1, FS -1, bottom +1, outer +1.');
    codeVersion = get_rankine_code_version();
    helperRelativeFile = fullfile('Source Code', '4.Potential', ...
        'get_rankine_source_orientation.m');
    helperFingerprintPresent = any(strcmp(codeVersion.relativeFiles, helperRelativeFile));
    rows(3) = make_after_result('T2_HELPER_IN_CODE_FINGERPRINT', ...
        double(~helperFingerprintPresent), 0.0, true, 'REPAIRED_PASS', ...
        'The production source-orientation helper participates in the solver code hash.');
    for testIndex = 1:numel(unaffectedIds)
        sourceRow = beforeTable(beforeTable.testId == unaffectedIds(testIndex), :);
        assert(height(sourceRow) == 1, 'CRESTU:AfterEvidenceSource', ...
            'Expected one fail-before row for %s.', unaffectedIds(testIndex));
        rows(3 + testIndex) = make_after_result(unaffectedIds(testIndex), ...
            sourceRow.relativeError, sourceRow.tolerance, false, ...
            'UNCHANGED_FAIL', ...
            'Unaffected analytic defect must remain currentPass=false after the T2-only repair.');
    end

    afterTable = struct2table(rows);
    csvFile = fullfile(evidenceDirectory, 'Orientation_Atomic_After_Results.csv');
    matFile = fullfile(evidenceDirectory, 'Orientation_Atomic_After_Results.mat');
    reportFile = fullfile(evidenceDirectory, 'ORIENTATION_ATOMIC_AFTER_EVIDENCE.md');
    writetable(afterTable, csvFile);
    metadata = struct('convention', orientation.convention, ...
        'componentNames', orientation.componentNames, ...
        'componentSigns', orientation.componentSigns, ...
        'orientationHash', orientation.signatureHash, ...
        'minimumWindingAlignment', orientation.minimumWindingAlignment, ...
        'codeFingerprint', codeVersion.fingerprint);
    save(matFile, 'afterTable', 'metadata', 'canonicalD', 'panels', ...
        'normals', 'centers');
    write_after_report(reportFile, afterTable, metadata);

    expectedMask = afterTable.currentPass == afterTable.expectedAfterPass;
    fprintf('[AUDIT] Orientation atomic after evidence: %d/%d rows match contract.\n', ...
        nnz(expectedMask), height(afterTable));
    fprintf('[AUDIT] Repaired mixed constant-field residual: %.16g\n', repairedResidual);
    fprintf('[AUDIT] Orientation hash: %s\n', orientation.signatureHash);
    fprintf('[AUDIT] CSV: %s\n', csvFile);
    fprintf('[AUDIT] MAT: %s\n', matFile);
    fprintf('[AUDIT] MD : %s\n', reportFile);
    assert(all(expectedMask), 'CRESTU:OrientationAtomicContract', ...
        'The atomic orientation after-state is broader or weaker than authorized.');
end

function row = empty_after_result()
    row = struct('testId', "", 'metric', NaN, 'tolerance', NaN, ...
        'currentPass', false, 'expectedAfterPass', false, ...
        'status', "", 'notes', "");
end

function row = make_after_result(testId, metric, tolerance, expectedPass, status, notes)
    row = empty_after_result();
    row.testId = string(testId);
    row.metric = metric;
    row.tolerance = tolerance;
    row.currentPass = metric <= tolerance;
    row.expectedAfterPass = expectedPass;
    row.status = string(status);
    row.notes = string(notes);
end

function panels = cube_panels_current_mixed(halfWidth)
    a = halfWidth;
    outward = zeros(6, 4, 3);
    outward(1, :, :) = [-a,-a,+a; +a,-a,+a; +a,+a,+a; -a,+a,+a];
    outward(2, :, :) = [-a,-a,-a; -a,+a,-a; +a,+a,-a; +a,-a,-a];
    outward(3, :, :) = [+a,-a,-a; +a,+a,-a; +a,+a,+a; +a,-a,+a];
    outward(4, :, :) = [-a,-a,-a; -a,-a,+a; -a,+a,+a; -a,+a,-a];
    outward(5, :, :) = [-a,+a,-a; -a,+a,+a; +a,+a,+a; +a,+a,-a];
    outward(6, :, :) = [-a,-a,-a; +a,-a,-a; +a,-a,+a; -a,-a,+a];
    panels = outward(:, [1,4,3,2], :);
    panels(1, :, :) = outward(1, :, :); % current FS stored/winding +ez
end

function [normals, centers] = panel_normals_and_centers(panels)
    panelCount = size(panels, 1);
    normals = zeros(panelCount, 3);
    centers = zeros(panelCount, 3);
    for panelIndex = 1:panelCount
        panel = reshape(panels(panelIndex, :, :), 4, 3);
        centers(panelIndex, :) = mean(panel, 1);
        areaVector = cross(panel(3, :) - panel(1, :), ...
            panel(4, :) - panel(2, :));
        normals(panelIndex, :) = areaVector / norm(areaVector);
    end
end

function D = assemble_canonical_double_layer(panels, normals, centers, signs)
    panelCount = size(panels, 1);
    D = zeros(panelCount);
    for collocationIndex = 1:panelCount
        for sourceIndex = 1:panelCount
            [~, dGdn] = rankine_panel_integrals( ...
                reshape(panels(sourceIndex, :, :), 4, 3), ...
                centers(collocationIndex, :), normals(sourceIndex, :));
            D(collocationIndex, sourceIndex) = ...
                -signs(sourceIndex) * dGdn / (4 * pi);
            if collocationIndex == sourceIndex
                D(collocationIndex, sourceIndex) = ...
                    D(collocationIndex, sourceIndex) + 0.5;
            end
        end
    end
end

function write_after_report(filename, afterTable, metadata)
    fileId = fopen(filename, 'w');
    if fileId < 0
        error('CRESTU:AfterReportWrite', 'Cannot write %s.', filename);
    end
    cleanup = onCleanup(@() fclose(fileId));
    fprintf(fileId, '# Double-Layer Orientation Atomic After Evidence\n\n');
    fprintf(fileId, '- Convention: `%s`\n', metadata.convention);
    fprintf(fileId, '- Component signs `[body, FS, bottom, outer]`: `%s`\n', ...
        mat2str(metadata.componentSigns));
    fprintf(fileId, '- Orientation hash: `%s`\n', metadata.orientationHash);
    fprintf(fileId, '- Minimum stored-normal/winding alignment: `%s`\n\n', ...
        mat2str(metadata.minimumWindingAlignment, 16));
    fprintf(fileId, '| Test | Metric | Tolerance | Current pass | Expected after pass | Status |\n');
    fprintf(fileId, '|---|---:|---:|---:|---:|---|\n');
    for rowIndex = 1:height(afterTable)
        fprintf(fileId, '| %s | %.16g | %.9g | %s | %s | %s |\n', ...
            afterTable.testId(rowIndex), afterTable.metric(rowIndex), ...
            afterTable.tolerance(rowIndex), string(afterTable.currentPass(rowIndex)), ...
            string(afterTable.expectedAfterPass(rowIndex)), afterTable.status(rowIndex));
    end
    fprintf(fileId, '\nOnly T2 is repaired. Rows marked `UNCHANGED_FAIL` remain physical failures, not successful tests.\n');
    clear cleanup;
end
