function result = P_FULL_M1_HANKEL_ABC_BENCHMARK()
%P_FULL_M1_HANKEL_ABC_BENCHMARK Full-normal-derivative Hankel ABC benchmark.
%
% This is deliberately self-contained.  It reads the frozen BMF directly,
% without read_bmf or any production geometry, mode, or ABC helper.  The
% three projections reported below are:
%   (1) full exact normal derivative, including (n_theta/r)*dtheta;
%   (2) exact axisymmetric geometric projection, omitting that term only;
%   (3) the current local production ABC projection.
%
% The theta-omission relative-error denominator is, panel by panel,
%   max(abs(fullExact), abs(k*Phi), eps).
% The baseline and total relative errors use analogous denominators based on
% their respective reference (axisymmetric and full) values.  The current
% production projection is retained as a known radial Hankel asymptotic
% closure; its difference from the exact axisymmetric projection is reported
% separately and is not treated as evidence that the closure is valid.

scriptDirectory = fileparts(mfilename('fullpath'));
if isempty(scriptDirectory)
    scriptDirectory = pwd;
end
bmfFile = fullfile(scriptDirectory, 'Frustum_Instrumentation', 'Cases', ...
    'Frustum_Fine_Top_Only_farfield.bmf');

[panelVertices, centers, normals, geometry] = local_read_frozen_bmf(bmfFile);
radius = hypot(centers(:, 1), centers(:, 2));
if any(~isfinite(radius)) || any(radius <= 0)
    error('P_FULL_M1_HANKEL_ABC:InvalidRadius', ...
        'All panel centers must have finite, positive cylindrical radius.');
end
er = [centers(:, 1) ./ radius, centers(:, 2) ./ radius, ...
    zeros(size(radius))];
etheta = [-centers(:, 2) ./ radius, centers(:, 1) ./ radius, ...
    zeros(size(radius))];
nr = sum(normals .* er, 2);
ntheta = sum(normals .* etheta, 2);
nz = normals(:, 3);

omega = 1.5;
g = 9.80665;
h = 50;
k = local_finite_depth_wavenumber(omega, g, h);
epsScale = eps;

% A: synthetic n_theta = 0.  This checks that the full and exact
% axisymmetric projections coincide to machine precision.  The production
% projection remains in the component report to expose its separate closure
% error.
modeA = local_mode_projection(centers, nr, zeros(size(ntheta)), nz, 1, ...
    k, h);
thetaZeroA = max(abs(modeA.fullExact - modeA.exactAxisymmetric));
rowA = local_summary('A', 'A_synthetic_ntheta0', 1, 0, modeA, ...
    local_components_finite(centers, normals, modeA), true, ...
    thetaZeroA <= 32 * epsScale * max(1, max(abs(modeA.fullExact))), ...
    ['Synthetic ntheta=0: full exact equals exact axisymmetric; ', ...
    'any production difference is only the known radial Hankel asymptotic closure.']);
% This is intentionally separate from passFlag: with ntheta=0, any
% remaining full-vs-production discrepancy is the pre-existing radial
% Hankel asymptotic closure error, not a theta omission.
rowA.productionReductionRelativeError = rowA.totalRelativeError;
rowA.productionReductionFlag = rowA.productionReductionRelativeError <= ...
    32 * epsScale;

% B: m = 0 on the actual geometry.  dtheta is identically zero, regardless
% of the actual n_theta values.
modeB = local_mode_projection(centers, nr, ntheta, nz, 0, k, h);
thetaZeroB = max(abs(modeB.fullExact - modeB.exactAxisymmetric));
rowB = local_summary('B', 'B_m0_actual_geometry', 0, ...
    max(abs(ntheta)), modeB, local_components_finite(centers, normals, modeB), ...
    all(nr < 0), thetaZeroB <= 32 * epsScale * max(1, max(abs(modeB.fullExact))), ...
    'm=0: theta contribution is exactly zero even with actual ntheta.');

% C: m = 1 on the actual geometry, followed by the requested synthetic
% ntheta magnitudes.  The synthetic orientation uses one actual panel as a
% representative: n_z is preserved, n_r is negative and rebuilt by the
% square-root closure so [n_r,n_theta,n_z] has unit norm.
modeCActual = local_mode_projection(centers, nr, ntheta, nz, 1, k, h);
actualFinite = all(isfinite([centers(:); normals(:); modeCActual.fullExact(:); ...
    modeCActual.exactAxisymmetric(:); modeCActual.currentProduction(:)]));
actualSign = all(nr < 0);
rowCActual = local_summary('C', 'C_m1_actual_geometry', 1, ...
    max(abs(ntheta)), modeCActual, actualFinite, actualSign, ...
    true, ...
    'm=1 actual geometry; winding normals and actual ntheta.');

representativeCandidates = find(isfinite(nr) & isfinite(nz) & (nr < 0));
if isempty(representativeCandidates)
    representativeIndex = 1;
else
    [~, localIndex] = max(abs(ntheta(representativeCandidates)));
    representativeIndex = representativeCandidates(localIndex);
end
nzRepresentative = nz(representativeIndex);
syntheticMagnitudes = [0, 0.002, 0.005, 0.01, ...
    0.01635907986834724, 0.02];
syntheticRows = repmat(rowCActual, numel(syntheticMagnitudes), 1);
syntheticComponents = cell(numel(syntheticMagnitudes), 1);
for q = 1:numel(syntheticMagnitudes)
    targetNtheta = syntheticMagnitudes(q);
    closureArgument = 1 - targetNtheta^2 - nzRepresentative^2;
    if closureArgument < 0
        nrSynthetic = NaN;
    else
        nrSynthetic = -sqrt(closureArgument);
    end
    % The representative center is used deliberately: this isolates the
    % n_theta magnitude experiment from any second geometric perturbation.
    representativeCenter = centers(representativeIndex, :);
    modeSynthetic = local_mode_projection(representativeCenter, nrSynthetic, ...
        targetNtheta, nzRepresentative, 1, k, h);
    syntheticFinite = all(isfinite([representativeCenter(:); nrSynthetic; ...
        targetNtheta; nzRepresentative; modeSynthetic.fullExact(:); ...
        modeSynthetic.exactAxisymmetric(:); modeSynthetic.currentProduction(:)]));
    syntheticNorm = nrSynthetic^2 + targetNtheta^2 + nzRepresentative^2;
    syntheticSign = isfinite(nrSynthetic) && nrSynthetic < 0 && ...
        abs(syntheticNorm - 1) <= 128 * epsScale;
    syntheticRows(q) = local_summary('C', ...
        sprintf('C_m1_synthetic_ntheta_%0.12g', targetNtheta), 1, ...
        targetNtheta, modeSynthetic, syntheticFinite, syntheticSign, true, ...
        'm=1 synthetic representative normal; n_z preserved, sqrt closure.');
    syntheticComponents{q} = modeSynthetic;
end

rows = [rowA; rowB; rowCActual; syntheticRows];
csvFile = fullfile(scriptDirectory, 'P_FULL_M1_HANKEL_ABC_BENCHMARK.csv');
local_write_csv(csvFile, rows);

result = struct();
result.BMFFile = bmfFile;
result.CSVFile = csvFile;
result.Parameters = struct('omega', omega, 'g', g, 'h', h, 'k', k, ...
    'timeConvention', 'exp(+i omega t)', ...
    'thetaRelativeDenominator', 'max(abs(fullExact),abs(k*Phi),eps)');
result.Geometry = geometry;
result.PanelVertices = panelVertices;
result.Centers = centers;
result.WindingNormals = normals;
result.Cylindrical = struct('radius', radius, 'nr', nr, ...
    'ntheta', ntheta, 'nz', nz, 'maxActualAbsNtheta', max(abs(ntheta)));
result.RepresentativeIndex = representativeIndex;
result.RepresentativeNz = nzRepresentative;
result.Rows = rows;
result.Components = struct('A', modeA, 'B', modeB, ...
    'CActual', modeCActual, 'CSynthetic', syntheticComponents);
result.Pass = all([rows.passFlag]);
end

function [vertices, centers, normals, geometry] = local_read_frozen_bmf(fileName)
% Read exactly the frozen BMF layout: four header lines, count, then 4x3
% coordinates per panel.  No production parser is called here.
fid = fopen(fileName, 'r');
if fid < 0
    error('P_FULL_M1_HANKEL_ABC:MissingBMF', 'Cannot open BMF: %s', fileName);
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
header = cell(4, 1);
for q = 1:4
    header{q} = fgetl(fid); %#ok<AGROW>
end
panelCount = sscanf(header{4}, '%d', 1);
if isempty(panelCount) || panelCount ~= 484
    error('P_FULL_M1_HANKEL_ABC:BMFPanelCount', ...
        'Expected 484 panels in frozen BMF, found %g.', panelCount);
end
vertices = zeros(panelCount, 4, 3);
for p = 1:panelCount
    for v = 1:4
        line = fgetl(fid);
        while ischar(line) && isempty(strtrim(line))
            line = fgetl(fid);
        end
        xyz = sscanf(line, '%f', 3);
        if numel(xyz) ~= 3
            error('P_FULL_M1_HANKEL_ABC:BMFVertex', ...
                'Panel %d vertex %d is not a 3-vector.', p, v);
        end
        vertices(p, v, :) = reshape(xyz, 1, 1, 3);
    end
end

centers = zeros(panelCount, 3);
normals = zeros(panelCount, 3);
for p = 1:panelCount
    panel = squeeze(vertices(p, :, :));
    uniquePanel = zeros(4, 3);
    uniqueCount = 0;
    for v = 1:4
        if uniqueCount == 0 || all(max(abs(uniquePanel(1:uniqueCount, :) - ...
                panel(v, :)), [], 2) > 1e-12)
            uniqueCount = uniqueCount + 1;
            uniquePanel(uniqueCount, :) = panel(v, :);
        end
    end
    uniquePanel = uniquePanel(1:uniqueCount, :);
    if uniqueCount < 3
        error('P_FULL_M1_HANKEL_ABC:DegeneratePanel', ...
            'Panel %d has fewer than three unique vertices.', p);
    end
    centers(p, :) = mean(uniquePanel, 1);
    found = false;
    for i = 1:(uniqueCount - 2)
        for j = (i + 1):(uniqueCount - 1)
            for ell = (j + 1):uniqueCount
                candidate = cross(uniquePanel(j, :) - uniquePanel(i, :), ...
                    uniquePanel(ell, :) - uniquePanel(i, :));
                candidateNorm = norm(candidate);
                if candidateNorm > 1e-14
                    normals(p, :) = candidate ./ candidateNorm;
                    found = true;
                    break;
                end
            end
            if found, break; end
        end
        if found, break; end
    end
    if ~found
        error('P_FULL_M1_HANKEL_ABC:DegeneratePanel', ...
            'Panel %d has no non-collinear winding triple.', p);
    end
end
geometry = struct('panelCount', panelCount, 'verticesPerPanel', 4, ...
    'uniqueVertexCounts', local_unique_counts(vertices));
end

function counts = local_unique_counts(vertices)
counts = zeros(size(vertices, 1), 1);
for p = 1:size(vertices, 1)
    panel = squeeze(vertices(p, :, :));
    uniquePanel = zeros(4, 3);
    n = 0;
    for v = 1:4
        if n == 0 || all(max(abs(uniquePanel(1:n, :) - panel(v, :)), [], 2) > 1e-12)
            n = n + 1;
            uniquePanel(n, :) = panel(v, :);
        end
    end
    counts(p) = n;
end
end

function k = local_finite_depth_wavenumber(omega, g, h)
target = omega^2 / g;
lo = 0;
hi = max(1, 2 * target);
while hi * tanh(hi * h) < target
    hi = 2 * hi;
end
for q = 1:100
    mid = 0.5 * (lo + hi);
    if mid * tanh(mid * h) < target
        lo = mid;
    else
        hi = mid;
    end
end
k = 0.5 * (lo + hi);
end

function components = local_mode_projection(centers, nr, ntheta, nz, m, k, h)
r = hypot(centers(:, 1), centers(:, 2));
theta = atan2(centers(:, 2), centers(:, 1));
z = centers(:, 3);
argument = k .* r;
hm = besselh(m, 2, argument);
hmNext = besselh(m + 1, 2, argument);
hmPrime = (m ./ argument) .* hm - hmNext;
phase = exp(1i * m .* theta);
zFactor = cosh(k .* (z + h)) ./ cosh(k .* h);
phi = hm .* phase .* zFactor;
dr = k .* hmPrime .* phase .* zFactor;
dtheta = 1i * m .* phi;
dz = k .* tanh(k .* (z + h)) .* phi;
fullExact = nr .* dr + (ntheta ./ r) .* dtheta + nz .* dz;
exactAxisymmetric = nr .* dr + nz .* dz;
currentProduction = ((-nr) .* (1i * k + 1 ./ (2 .* r)) + ...
    nz .* k .* tanh(k .* (z + h))) .* phi;

% The theta-only denominator is the frozen benchmark convention.
thetaDenominator = max(max(abs(fullExact), abs(k .* phi)), eps);
thetaAbsolute = abs(fullExact - exactAxisymmetric);
thetaRelative = thetaAbsolute ./ thetaDenominator;
baselineDenominator = max(max(abs(exactAxisymmetric), abs(k .* phi)), eps);
baselineAbsolute = abs(exactAxisymmetric - currentProduction);
baselineRelative = baselineAbsolute ./ baselineDenominator;
totalDenominator = max(max(abs(fullExact), abs(k .* phi)), eps);
totalAbsolute = abs(fullExact - currentProduction);
totalRelative = totalAbsolute ./ totalDenominator;
components = struct('m', m, 'Phi', phi, 'dr', dr, 'dtheta', dtheta, 'dz', dz, ...
    'fullExact', fullExact, 'exactAxisymmetric', exactAxisymmetric, ...
    'currentProduction', currentProduction, 'thetaAbsolute', thetaAbsolute, ...
    'thetaRelative', thetaRelative, 'baselineAbsolute', baselineAbsolute, ...
    'baselineRelative', baselineRelative, 'totalAbsolute', totalAbsolute, ...
    'totalRelative', totalRelative);
end

function flag = local_components_finite(centers, normals, components)
flag = all(isfinite([centers(:); normals(:); components.Phi(:); ...
    components.dr(:); components.dtheta(:); components.dz(:); ...
    components.fullExact(:); components.exactAxisymmetric(:); ...
    components.currentProduction(:); components.thetaAbsolute(:); ...
    components.thetaRelative(:); components.baselineRelative(:); ...
    components.totalRelative(:)]));
end

function row = local_summary(section, caseName, m, nthetaValue, components, ...
    finiteFlag, signFlag, identityFlag, description)
thetaAbs = abs(components.thetaAbsolute(:));
thetaRel = abs(components.thetaRelative(:));
baseRel = abs(components.baselineRelative(:));
totalRel = abs(components.totalRelative(:));
row = struct();
row.section = section;
row.caseName = caseName;
row.m = m;
row.ntheta = nthetaValue;
row.maxThetaAbsolute = max(thetaAbs);
row.medianThetaAbsolute = median(thetaAbs);
row.p95ThetaAbsolute = local_percentile(thetaAbs, 95);
row.maxThetaRelative = max(thetaRel);
row.medianThetaRelative = median(thetaRel);
row.p95ThetaRelative = local_percentile(thetaRel, 95);
row.baselineAsymptoticRelativeError = max(baseRel);
row.medianBaselineAsymptoticRelativeError = median(baseRel);
row.p95BaselineAsymptoticRelativeError = local_percentile(baseRel, 95);
row.totalRelativeError = max(totalRel);
row.medianTotalRelativeError = median(totalRel);
row.p95TotalRelativeError = local_percentile(totalRel, 95);
row.finiteFlag = logical(finiteFlag);
row.signFlag = logical(signFlag);
row.identityFlag = logical(identityFlag);
row.productionReductionFlag = false;
row.productionReductionRelativeError = NaN;
row.passFlag = logical(finiteFlag && signFlag && identityFlag);
row.description = description;
end

function value = local_percentile(values, percent)
values = sort(values(isfinite(values)));
if isempty(values)
    value = NaN;
    return;
end
if numel(values) == 1
    value = values(1);
    return;
end
position = 1 + (numel(values) - 1) * percent / 100;
lower = floor(position);
upper = ceil(position);
if lower == upper
    value = values(lower);
else
    value = values(lower) + (position - lower) * (values(upper) - values(lower));
end
end

function local_write_csv(fileName, rows)
fid = fopen(fileName, 'w');
if fid < 0
    error('P_FULL_M1_HANKEL_ABC:CSV', 'Cannot write CSV: %s', fileName);
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, ['section,case,m,ntheta,max_theta_absolute,median_theta_absolute,', ...
    'p95_theta_absolute,max_theta_relative,median_theta_relative,p95_theta_relative,', ...
    'baseline_asymptotic_relative_error,median_baseline_asymptotic_relative_error,', ...
    'p95_baseline_asymptotic_relative_error,total_relative_error,', ...
    'median_total_relative_error,p95_total_relative_error,finite_flag,sign_flag,', ...
    'identity_flag,production_reduction_flag,production_reduction_relative_error,', ...
    'pass_flag,description\n']);
for q = 1:numel(rows)
    description = strrep(rows(q).description, '"', '""');
    fprintf(fid, ['%s,"%s",%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,', ...
        '%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%d,%d,%d,%d,%.17g,%d,"%s"\n'], ...
        rows(q).section, rows(q).caseName, rows(q).m, rows(q).ntheta, ...
        rows(q).maxThetaAbsolute, rows(q).medianThetaAbsolute, ...
        rows(q).p95ThetaAbsolute, rows(q).maxThetaRelative, ...
        rows(q).medianThetaRelative, rows(q).p95ThetaRelative, ...
        rows(q).baselineAsymptoticRelativeError, ...
        rows(q).medianBaselineAsymptoticRelativeError, ...
        rows(q).p95BaselineAsymptoticRelativeError, rows(q).totalRelativeError, ...
        rows(q).medianTotalRelativeError, rows(q).p95TotalRelativeError, ...
        rows(q).finiteFlag, rows(q).signFlag, rows(q).identityFlag, ...
        rows(q).productionReductionFlag, rows(q).productionReductionRelativeError, ...
        rows(q).passFlag, description);
end
end
