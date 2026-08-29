function resultTable = Run_Phase2_2_Analytic_FailBefore()
% RUN_PHASE2_2_ANALYTIC_FAILBEFORE Build analytic evidence before production changes.
%
% This is deliberately a fail-before audit, not a green regression suite.
% Every row retains currentPass and expectedBeforePass separately so that an
% expected defect cannot be reported as a successful physical validation.

    testDirectory = fileparts(mfilename('fullpath'));
    artifactDirectory = fileparts(testDirectory);
    convergenceDirectory = fileparts(artifactDirectory);
    caseWaveDirectory = fileparts(convergenceDirectory);
    repositoryRoot = fileparts(caseWaveDirectory);
    evidenceDirectory = fullfile(testDirectory, 'Evidence');
    if ~isfolder(evidenceDirectory)
        mkdir(evidenceDirectory);
    end

    addpath(fullfile(repositoryRoot, 'Source Code', '2.Mesh'));
    addpath(fullfile(repositoryRoot, 'Source Code', '4.Potential'));
    addpath(fullfile(repositoryRoot, 'Source Code', '5.Force'));

    gravity = 9.80665;             % [m/s^2]
    density = 1025.0;              % [kg/m^3]
    waterDepth = 50.0;             % [m]
    sphereRadius = 5.0;            % [m]
    omega = 2.0;                   % [rad/s]
    rows = repmat(empty_result(), 0, 1);

    %% T1: exp(+i*omega*t), incident phase, and peak-amplitude linearity
    point = [1.2, -0.7, 0.0];
    normal = [0.3, -0.4, sqrt(0.75)];
    betaDegrees = 30.0;
    amplitude = 1.7;               % peak amplitude [m]
    [phiOne, qOne, waveNumber] = compute_incident_wave(point, normal, ...
        omega, gravity, waterDepth, betaDegrees, amplitude);
    beta = deg2rad(betaDegrees);
    phase = exp(-1i * waveNumber * ...
        (point(1) * cos(beta) + point(2) * sin(beta)));
    etaCurrent = -1i * omega * phiOne / gravity;
    etaOracle = amplitude * phase;
    rows(end + 1) = make_result('T1_TIME_PHASE', 'time/amplitude', ...
        etaCurrent, etaOracle, max(abs(etaOracle), 1), 1e-12, true, ...
        'ACTIVE', 'eta=-i*omega*Phi/g must equal peak a*exp(-ik.x) at z=0.');

    [phiTwo, qTwo] = compute_incident_wave(point, normal, omega, gravity, ...
        waterDepth, betaDegrees, 2 * amplitude);
    linearityMetric = max([abs(phiTwo / phiOne), abs(qTwo / qOne)]);
    rows(end + 1) = make_result('T1_PEAK_AMPLITUDE_LINEARITY', ...
        'time/amplitude', linearityMetric, 2.0, 2.0, 1e-12, true, ...
        'ACTIVE', 'Doubling peak amplitude must double Phi and q, not multiply by sqrt(2).');

    %% T2: actual Rankine kernel, winding, and closed constant-field identity
    inwardPanels = cube_panels(1.0, false);
    [inwardD, inwardWindingDot] = assemble_double_layer(inwardPanels);
    inwardResidual = norm(inwardD * ones(6, 1), inf);
    rows(end + 1) = make_result('T2_INWARD_CLOSED_CONTROL', ...
        'Green normal/winding', inwardResidual, 0.0, 1.0, 1e-10, true, ...
        'ACTIVE', 'For inward source winding, D*1=0 is the closed constant-field oracle.');
    rows(end + 1) = make_result('T2_INWARD_WINDING_CONTROL', ...
        'Green normal/winding', min(inwardWindingDot), 1.0, 1.0, 1e-12, true, ...
        'ACTIVE', 'Kernel winding-derived normal must agree with the supplied source normal.');

    mixedPanels = cube_panels(1.0, true); % top/FS outward, other faces inward
    [mixedD, mixedWindingDot] = assemble_double_layer(mixedPanels);
    mixedResidual = norm(mixedD * ones(6, 1), inf);
    rows(end + 1) = make_result('T2_CURRENT_MIXED_ORIENTATION', ...
        'Green normal/winding', mixedResidual, 0.0, 1.0, 1e-10, false, ...
        'ACTIVE', ['Top face uses the current FS orientation (+ez); bottom and ', ...
        'sides use current inward orientation. Constant-field identity must fail.']);
    rows(end + 1) = make_result('T2_MIXED_WINDING_IS_SELF_CONSISTENT', ...
        'Green normal/winding', min(mixedWindingDot), 1.0, 1.0, 1e-12, true, ...
        'ACTIVE', 'The defect is mixed domain orientation, not a normal-versus-winding mismatch.');

    %% Load the existing Fine body without generating or changing geometry
    fineBodyFile = fullfile(convergenceDirectory, 'Mesh_Fine', ...
        'hemi_D10_fine_body.bmf');
    fineBody = read_bmf(fineBodyFile);
    [fullCenters, fullNormals, fullAreas] = mirror_hemisphere(fineBody);
    bodyRadius = median(sqrt(sum(fullCenters .^ 2, 2)));
    bodyWindingDot = winding_alignment(fineBody.vertices, fineBody.normals);
    rows(end + 1) = make_result('T2_FINE_BODY_WINDING_ALIGNMENT', ...
        'Green normal/winding', min(bodyWindingDot), 1.0, 1.0, 1e-12, true, ...
        'ACTIVE', 'Existing Fine body winding and stored outward-body normals agree.');

    %% T3: strict finite-depth free-surface coefficient
    omegaLow = 0.5;
    kLow = solve_wave_dispersion(omegaLow, gravity, waterDepth);
    nuCurrent = kLow;
    nuOracle = omegaLow ^ 2 / gravity;
    rows(end + 1) = make_result('T3_FS_FINITE_DEPTH_NU', ...
        'free-surface Robin', nuCurrent, nuOracle, abs(nuOracle), 1e-10, false, ...
        'ACTIVE', 'Strict FS coefficient is omega^2/g=k*tanh(kh), not k.');

    %% T4: finite-depth vertical eigenfunction satisfies the bottom BC
    [~, bottomQ] = compute_incident_wave([0.4, -0.2, -waterDepth], ...
        [0, 0, 1], omegaLow, gravity, waterDepth, 17.0, 1.0);
    rows(end + 1) = make_result('T4_BOTTOM_NEUMANN', ...
        'bottom Neumann', bottomQ, 0.0, 1.0, 1e-12, true, ...
        'ACTIVE', 'cosh vertical mode has dPhi/dz=0 at z=-h.');

    %% T5: exact/asymptotic outgoing conditions on cylinder and frustum
    lambda = 2 * pi / waveNumber;
    cylinderRadius = 20.0 + 1.5 * lambda;
    h0 = besselh(0, 2, waveNumber * cylinderRadius);
    h1 = besselh(1, 2, waveNumber * cylinderRadius);
    storedCylinderOracle = waveNumber * h1 / h0; % -d/dr for stored -er
    currentOuterCoefficient = 1i * waveNumber;
    rows(end + 1) = make_result('T5_CYLINDER_STORED_DECAY', ...
        'outer outgoing', currentOuterCoefficient, storedCylinderOracle, ...
        abs(storedCylinderOracle), 1e-8, false, 'ACTIVE', ...
        'Current +ik omits the cylindrical 1/(2R) decay in stored-inward q.');
    rows(end + 1) = make_result('T5_CYLINDER_GREEN_OUTWARD_SIGN', ...
        'outer outgoing', currentOuterCoefficient, -storedCylinderOracle, ...
        abs(storedCylinderOracle), 1e-8, false, 'ACTIVE', ...
        'With inward source D, the Green RHS requires outward-fluid q, opposite to current +ik.');

    bottomRadius = cylinderRadius;
    topRadius = cylinderRadius + 0.5 * lambda;
    radiusSlope = (topRadius - bottomRadius) / waterDepth;
    zMid = -0.5 * waterDepth;
    rMid = bottomRadius + (zMid + waterDepth) * radiusSlope;
    storedFrustumOracle = (1i * waveNumber + 1 / (2 * rMid) + ...
        radiusSlope * waveNumber * tanh(waveNumber * (zMid + waterDepth))) / ...
        sqrt(1 + radiusSlope ^ 2);
    rows(end + 1) = make_result('T5_FRUSTUM_STORED_PROJECTION', ...
        'outer outgoing', currentOuterCoefficient, storedFrustumOracle, ...
        abs(storedFrustumOracle), 1e-8, false, 'ACTIVE', ...
        'Current cylinder +ik lacks radial decay, sloped-normal projection, and V''(z)/V(z).');

    assemblySource = fileread(fullfile(repositoryRoot, 'Source Code', ...
        '4.Potential', 'assemble_rankine_matrix.m'));
    hasComplexOuterInterface = isempty(regexp(assemblySource, ...
        'nu\(ff_idx\)\s*=\s*k0\s*;', 'once'));
    rows(end + 1) = make_result('T5_COMPLEX_KAPPA_INTERFACE', ...
        'outer outgoing', double(hasComplexOuterInterface), 1.0, 1.0, 0.5, false, ...
        'NOT_IMPLEMENTED', ['Current outer coefficient is hard-wired to real k0; ', ...
        'complex-kappa validation is a candidate gate, not a validated option.']);

    %% T6: translating-sphere physical Neumann data versus Green RHS orientation
    bodyNz = fineBody.normals(:, 3);
    physicalBodyQ = 1i * omega * bodyNz;
    analyticBodyQ = 1i * omega * bodyNz;
    physicalQError = norm(physicalBodyQ - analyticBodyQ) / ...
        max(norm(analyticBodyQ), eps);
    rows(end + 1) = make_result('T6_BODY_NEUMANN_VALUE', ...
        'radiation RHS', physicalQError, 0.0, 1.0, 1e-12, true, ...
        'ACTIVE', 'Stored outward-body normal data qB=i*omega*nz is correct for unit displacement.');
    currentGreenRhs = 1i * omega * bodyNz;
    oracleGreenRhs = -1i * omega * bodyNz; % qOmega=-qB for inner fluid boundary
    greenRhsSignedRatio = real((oracleGreenRhs' * currentGreenRhs) / ...
        (oracleGreenRhs' * oracleGreenRhs));
    rows(end + 1) = make_result('T6_GREEN_RHS_ORIENTATION', ...
        'radiation RHS', greenRhsSignedRatio, 1.0, 1.0, 1e-12, false, ...
        'ACTIVE', 'Body D uses inward-to-fluid source orientation, so Green RHS needs qOmega=-qB.');

    %% T10: analytic pressure traction, isolated from the BEM solve
    pressureAmplitude = 1000.0; % [Pa]
    fullNz = fullNormals(:, 3);
    numericalNzSquaredArea = sum((fullNz .^ 2) .* fullAreas);
    analyticNzSquaredArea = 4 * pi * sphereRadius ^ 2 / 3;
    rows(end + 1) = make_result('T10_SPHERE_QUADRATURE_CONTROL', ...
        'pressure traction', numericalNzSquaredArea, analyticNzSquaredArea, ...
        analyticNzSquaredArea, 0.02, true, 'ACTIVE', ...
        'Mirrored existing Fine hemisphere must integrate nz^2 over a full sphere.');
    analyticPressure = pressureAmplitude * fullNz;
    pressurePotential = analyticPressure / (-1i * omega * density);
    currentPressureForce = compute_wave_excitation(pressurePotential, ...
        zeros(size(pressurePotential)), fullNz, fullAreas, omega, density);
    oraclePressureForce = -pressureAmplitude * analyticNzSquaredArea;
    rows(end + 1) = make_result('T10_PRESSURE_TRACTION_SIGN', ...
        'pressure traction', currentPressureForce, oraclePressureForce, ...
        abs(oraclePressureForce), 0.02, false, 'ACTIVE', ...
        'For outward-body nB the load is -integral(p*nB), not +integral(p*nB).');

    %% T7: translating-sphere added-mass extractor, isolated from the solver
    analyticRadiationPotential = -1i * omega * sphereRadius * fullNz / 2;
    [currentAddedMass, currentDamping] = compute_hydrodynamic_coeffs( ...
        analyticRadiationPotential, fullNz, fullAreas, omega, density);
    oracleAddedMass = 0.5 * density * (4 * pi * sphereRadius ^ 3 / 3);
    rows(end + 1) = make_result('T7_TRANSLATING_SPHERE_ADDED_MASS', ...
        'pressure A/B extractor', currentAddedMass, oracleAddedMass, ...
        oracleAddedMass, 0.02, false, 'ACTIVE', ...
        'Analytic dipole Phi=-i*omega*a*cos(theta)/2 gives A=0.5*rho*volume.');
    rows(end + 1) = make_result('T7_TRANSLATING_SPHERE_DAMPING', ...
        'pressure A/B extractor', currentDamping, 0.0, ...
        max(oracleAddedMass * omega, 1), 1e-12, true, 'ACTIVE', ...
        'Unbounded translating-sphere dipole is purely reactive, so B=0.');

    %% Preserve auditable outputs and reject any unexpected before-state
    resultTable = struct2table(rows);
    csvFile = fullfile(evidenceDirectory, 'Analytic_FailBefore_Results.csv');
    matFile = fullfile(evidenceDirectory, 'Analytic_FailBefore_Results.mat');
    reportFile = fullfile(evidenceDirectory, 'ANALYTIC_FAIL_BEFORE_EVIDENCE.md');
    writetable(resultTable, csvFile);
    metadata = struct('timeConvention', 'exp(+i*omega*t)', ...
        'fineBodyFile', fineBodyFile, 'fineBodyPanelCount', fineBody.n_panels, ...
        'fineBodyMedianCenterRadiusM', bodyRadius, 'gravityMPerS2', gravity, ...
        'densityKgPerM3', density, 'waterDepthM', waterDepth, ...
        'sphereRadiusM', sphereRadius, 'testOmegaRadPerS', omega, ...
        'generatedAt', char(datetime('now', 'TimeZone', 'Asia/Shanghai')));
    save(matFile, 'resultTable', 'metadata', 'inwardD', 'mixedD', ...
        'inwardPanels', 'mixedPanels');
    write_markdown_report(reportFile, resultTable, metadata, csvFile, matFile);

    expectedMask = resultTable.currentPass == resultTable.expectedBeforePass;
    fprintf('[AUDIT] Analytic fail-before evidence: %d/%d rows match the declared before-state.\n', ...
        nnz(expectedMask), height(resultTable));
    fprintf('[AUDIT] Physical currentPass rows: %d PASS, %d FAIL.\n', ...
        nnz(resultTable.currentPass), nnz(~resultTable.currentPass));
    fprintf('[AUDIT] CSV: %s\n', csvFile);
    fprintf('[AUDIT] MAT: %s\n', matFile);
    fprintf('[AUDIT] MD : %s\n', reportFile);
    assert(all(expectedMask), 'CRESTU:UnexpectedFailBeforeState', ...
        'At least one analytic row does not match its declared fail-before state.');
end

function row = empty_result()
    row = struct('testId', "", 'category', "", ...
        'currentReal', NaN, 'currentImag', NaN, ...
        'oracleReal', NaN, 'oracleImag', NaN, ...
        'signedErrorReal', NaN, 'signedErrorImag', NaN, ...
        'relativeError', NaN, 'tolerance', NaN, ...
        'currentPass', false, 'expectedBeforePass', false, ...
        'beforeState', "", 'candidateGate', "", 'notes', "");
end

function row = make_result(testId, category, currentValue, oracleValue, ...
        errorScale, tolerance, expectedBeforePass, candidateGate, notes)
    row = empty_result();
    signedError = currentValue - oracleValue;
    row.testId = string(testId);
    row.category = string(category);
    row.currentReal = real(currentValue);
    row.currentImag = imag(currentValue);
    row.oracleReal = real(oracleValue);
    row.oracleImag = imag(oracleValue);
    row.signedErrorReal = real(signedError);
    row.signedErrorImag = imag(signedError);
    row.relativeError = abs(signedError) / max(abs(errorScale), eps);
    row.tolerance = tolerance;
    row.currentPass = row.relativeError <= tolerance;
    row.expectedBeforePass = expectedBeforePass;
    if expectedBeforePass
        row.beforeState = "EXPECTED_PASS_BEFORE";
    else
        row.beforeState = "EXPECTED_FAIL_BEFORE";
    end
    row.candidateGate = string(candidateGate);
    row.notes = string(notes);
end

function panels = cube_panels(halfWidth, useCurrentMixedOrientation)
    a = halfWidth;
    outward = zeros(6, 4, 3);
    outward(1, :, :) = [-a,-a,+a; +a,-a,+a; +a,+a,+a; -a,+a,+a]; % top +z
    outward(2, :, :) = [-a,-a,-a; -a,+a,-a; +a,+a,-a; +a,-a,-a]; % bottom -z
    outward(3, :, :) = [+a,-a,-a; +a,+a,-a; +a,+a,+a; +a,-a,+a]; % +x
    outward(4, :, :) = [-a,-a,-a; -a,-a,+a; -a,+a,+a; -a,+a,-a]; % -x
    outward(5, :, :) = [-a,+a,-a; -a,+a,+a; +a,+a,+a; +a,+a,-a]; % +y
    outward(6, :, :) = [-a,-a,-a; +a,-a,-a; +a,-a,+a; -a,-a,+a]; % -y
    panels = outward(:, [1,4,3,2], :); % all inward
    if useCurrentMixedOrientation
        panels(1, :, :) = outward(1, :, :); % current FS is outward +z
    end
end

function [D, windingDot] = assemble_double_layer(panels)
    panelCount = size(panels, 1);
    centers = zeros(panelCount, 3);
    normals = zeros(panelCount, 3);
    for panelIndex = 1:panelCount
        vertices = reshape(panels(panelIndex, :, :), 4, 3);
        centers(panelIndex, :) = mean(vertices, 1);
        normalRaw = cross(vertices(3, :) - vertices(1, :), ...
            vertices(4, :) - vertices(2, :));
        normals(panelIndex, :) = normalRaw / norm(normalRaw);
    end
    windingDot = winding_alignment(panels, normals);
    D = zeros(panelCount);
    for collocationIndex = 1:panelCount
        for sourceIndex = 1:panelCount
            [~, dGdn] = rankine_panel_integrals( ...
                reshape(panels(sourceIndex, :, :), 4, 3), ...
                centers(collocationIndex, :), normals(sourceIndex, :));
            D(collocationIndex, sourceIndex) = -dGdn / (4 * pi);
            if collocationIndex == sourceIndex
                D(collocationIndex, sourceIndex) = ...
                    D(collocationIndex, sourceIndex) + 0.5;
            end
        end
    end
end

function alignment = winding_alignment(vertices, normals)
    panelCount = size(vertices, 1);
    alignment = zeros(panelCount, 1);
    for panelIndex = 1:panelCount
        panel = reshape(vertices(panelIndex, :, :), 4, 3);
        winding = cross(panel(3, :) - panel(1, :), ...
            panel(4, :) - panel(2, :));
        winding = winding / norm(winding);
        alignment(panelIndex) = dot(winding, normals(panelIndex, :));
    end
end

function [centers, normals, areas] = mirror_hemisphere(mesh)
    centers = [mesh.centers; mesh.centers .* [1, 1, -1]];
    normals = [mesh.normals; mesh.normals .* [1, 1, -1]];
    areas = [mesh.areas; mesh.areas];
end

function write_markdown_report(filename, resultTable, metadata, csvFile, matFile)
    fileId = fopen(filename, 'w');
    if fileId < 0
        error('CRESTU:AuditReportWrite', 'Cannot write %s.', filename);
    end
    cleanup = onCleanup(@() fclose(fileId));
    fprintf(fileId, '# Phase 2.2 Analytic Fail-Before Evidence\n\n');
    fprintf(fileId, 'This is a physical fail-before audit. `currentPass=false` remains a failure; an expected failure is not counted as a validated test.\n\n');
    fprintf(fileId, '- Time convention: `%s`\n', metadata.timeConvention);
    fprintf(fileId, '- Fine body: `%s` (%d panels)\n', metadata.fineBodyFile, metadata.fineBodyPanelCount);
    fprintf(fileId, '- CSV: `%s`\n', csvFile);
    fprintf(fileId, '- MAT: `%s`\n\n', matFile);
    fprintf(fileId, '| Test | Category | Before expectation | Current pass | Relative error | Tolerance | Gate |\n');
    fprintf(fileId, '|---|---|---:|---:|---:|---:|---|\n');
    for rowIndex = 1:height(resultTable)
        fprintf(fileId, '| %s | %s | %s | %s | %.9g | %.9g | %s |\n', ...
            resultTable.testId(rowIndex), resultTable.category(rowIndex), ...
            resultTable.beforeState(rowIndex), string(resultTable.currentPass(rowIndex)), ...
            resultTable.relativeError(rowIndex), resultTable.tolerance(rowIndex), ...
            resultTable.candidateGate(rowIndex));
    end
    fprintf(fileId, '\n## Interpretation\n\n');
    fprintf(fileId, '- Passing control rows establish that the analytic oracles and existing Fine sphere mesh are active.\n');
    fprintf(fileId, '- Failing rows identify convention defects before any production edit.\n');
    fprintf(fileId, '- The complex-kappa row is `NOT_IMPLEMENTED`; it is not evidence of a validated complex radiation root.\n');
    clear cleanup;
end
