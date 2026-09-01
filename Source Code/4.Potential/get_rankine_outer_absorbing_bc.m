function outerBC = get_rankine_outer_absorbing_bc(centers, normals, stats, ...
        omega, gravity, waterDepth, cfg, physicalK, nominalMu)
% GET_RANKINE_OUTER_ABSORBING_BC Build a first-order local asymptotic ABC.
% This is not an exact Dirichlet-to-Neumann operator.

    policy = outer_policy();
    if nargin == 0
        outerBC = struct('enabled', false, 'policy', policy, ...
            'gamma', complex(zeros(0, 1)), 'cacheSpecification', policy);
        return;
    end
    validateattributes(omega, {'numeric'}, {'scalar','real','positive','finite'});
    validateattributes(gravity, {'numeric'}, {'scalar','real','positive','finite'});
    validateattributes(waterDepth, {'numeric'}, {'scalar','real','positive','finite'});
    validateattributes(physicalK, {'numeric'}, {'scalar','real','positive','finite'});
    validateattributes(nominalMu, {'numeric'}, {'scalar','real','nonnegative','finite'});
    centers = reshape(centers, [], 3);
    normals = reshape(normals, [], 3);
    assert(size(centers, 1) == size(normals, 1), ...
        'CRESTU:OuterABCGeometry', 'Centers and normals have inconsistent rows.');

    bodyCount = stats.total_body_panels;
    freeSurfaceCount = stats.fs_panels;
    bottomCount = stats.seabed_panels;
    outerCount = stats.farfield_panels;
    expectedCount = bodyCount + freeSurfaceCount + bottomCount + outerCount;
    assert(size(centers, 1) == expectedCount, ...
        'CRESTU:OuterABCCounts', 'Boundary counts do not match the merged geometry.');
    if outerCount == 0
        cacheSpecification = policy;
        cacheSpecification.enabled = false;
        outerBC = struct('enabled', false, 'policy', policy, ...
            'gamma', complex(zeros(0, 1)), ...
            'cacheSpecification', cacheSpecification);
        return;
    end
    assert(freeSurfaceCount > 0, 'CRESTU:OuterABCFreeSurfaceRing', ...
        'An outer absorbing boundary requires a resolved free-surface outer ring.');

    freeSurfaceRows = bodyCount + (1:freeSurfaceCount);
    outerRows = bodyCount + freeSurfaceCount + bottomCount + (1:outerCount);
    freeSurfaceRadii = hypot(centers(freeSurfaceRows, 1), centers(freeSurfaceRows, 2));
    actualRingRadius = max(freeSurfaceRadii);
    spongeStart = cfg.fs.r_inner;
    spongeWidth = cfg.fs.r_outer - cfg.fs.r_inner;
    if isfield(cfg, 'phase2_2_controls') && ...
            isfield(cfg.phase2_2_controls, 'effective')
        spongeStart = cfg.phase2_2_controls.effective.spongeStartRadiusM;
        spongeWidth = cfg.phase2_2_controls.effective.spongeWidthM;
    end
    assert(isfinite(spongeWidth) && spongeWidth > 0, ...
        'CRESTU:OuterABCSpongeWidth', 'Sponge width must be positive and finite.');
    normalizedRingPosition = max(0, min(1, ...
        (actualRingRadius - spongeStart) / spongeWidth));
    actualMu = nominalMu * normalizedRingPosition ^ 2;

    nu0 = omega ^ 2 / gravity;
    [kappa, continuation] = continue_complex_kappa( ...
        physicalK, waterDepth, nu0, actualMu, policy);
    outerCenters = centers(outerRows, :);
    outerNormals = normals(outerRows, :);
    radius = hypot(outerCenters(:, 1), outerCenters(:, 2));
    assert(all(isfinite(radius) & radius > 0), ...
        'CRESTU:OuterABCRadius', 'Every outer source panel must have r>0.');
    radialNormal = (outerNormals(:, 1) .* outerCenters(:, 1) + ...
        outerNormals(:, 2) .* outerCenters(:, 2)) ./ radius;
    tangentialNormal = (-outerNormals(:, 1) .* outerCenters(:, 2) + ...
        outerNormals(:, 2) .* outerCenters(:, 1)) ./ radius;
    verticalNormal = outerNormals(:, 3);
    assert(all(radialNormal < 0), 'CRESTU:OuterABCRadialNormal', ...
        'Stored outer normals must have strictly inward radial projection.');
    assert(max(abs(tangentialNormal)) <= policy.normalTangentialTolerance, ...
        'CRESTU:OuterABCTangentialNormal', ...
        'Outer normals violate the axisymmetric local-ABC assumption.');

    radialCoefficient = 1i * kappa + 1 ./ (2 * radius);
    verticalCoefficient = kappa * tanh(kappa * (outerCenters(:, 3) + waterDepth));
    gamma = (-radialNormal) .* radialCoefficient + ...
        verticalNormal .* verticalCoefficient;
    assert(all(isfinite(real(gamma)) & isfinite(imag(gamma))), ...
        'CRESTU:OuterABCGammaFinite', 'Outer ABC gamma contains NaN or Inf.');

    modeErrors = zeros(2, 3);
    for mode = 0:1
        exactGamma = exact_cylinder_gamma(mode, kappa, radius);
        relativeError = abs(radialCoefficient - exactGamma) ./ ...
            max(abs(exactGamma), eps);
        modeErrors(mode + 1, :) = [min(relativeError), ...
            median(relativeError), max(relativeError)];
    end
    gammaHash = numeric_sha256([real(gamma); imag(gamma); radius; ...
        radialNormal; tangentialNormal; verticalNormal]);
    cacheSpecification = policy;
    cacheSpecification.enabled = true;
    cacheSpecification.physicalKPerM = physicalK;
    cacheSpecification.kappaRealPerM = real(kappa);
    cacheSpecification.kappaImagPerM = imag(kappa);
    cacheSpecification.nu0PerM = nu0;
    cacheSpecification.nominalMu = nominalMu;
    cacheSpecification.spongeStartRadiusM = spongeStart;
    cacheSpecification.spongeWidthM = spongeWidth;
    cacheSpecification.spongeEndRadiusM = spongeStart + spongeWidth;
    cacheSpecification.actualOuterRingMu = actualMu;
    cacheSpecification.actualOuterFreeSurfaceRingRadiusM = actualRingRadius;
    cacheSpecification.normalizedOuterRingPosition = normalizedRingPosition;
    cacheSpecification.continuationAcceptedSteps = continuation.acceptedSteps;
    cacheSpecification.continuationRefinementLevel = continuation.refinementLevel;
    cacheSpecification.newtonTotalIterations = continuation.totalIterations;
    cacheSpecification.newtonMaximumIterations = continuation.maximumIterations;
    cacheSpecification.dispersionRelativeResidual = continuation.relativeResidual;
    cacheSpecification.minimumContinuationRealKappa = continuation.minimumRealPart;
    cacheSpecification.maximumContinuationImagKappa = continuation.maximumImagPart;
    cacheSpecification.radiusRangeM = [min(radius), max(radius)];
    cacheSpecification.radialNormalRange = [min(radialNormal), max(radialNormal)];
    cacheSpecification.verticalNormalRange = [min(verticalNormal), max(verticalNormal)];
    cacheSpecification.maximumAbsTangentialNormal = max(abs(tangentialNormal));
    cacheSpecification.physicalKRadiusRange = physicalK * [min(radius), max(radius)];
    cacheSpecification.absKappaRadiusRange = abs(kappa) * [min(radius), max(radius)];
    cacheSpecification.m0AsymptoticErrorMinMedianMax = modeErrors(1, :);
    cacheSpecification.m1AsymptoticErrorMinMedianMax = modeErrors(2, :);
    cacheSpecification.gammaRealRangePerM = [min(real(gamma)), max(real(gamma))];
    cacheSpecification.gammaImagRangePerM = [min(imag(gamma)), max(imag(gamma))];
    cacheSpecification.gammaAbsRangePerM = [min(abs(gamma)), max(abs(gamma))];
    cacheSpecification.gammaSHA256 = gammaHash;
    outerBC = struct('enabled', true, 'policy', policy, ...
        'gamma', gamma, 'kappa', kappa, 'physicalK', physicalK, ...
        'nominalMu', nominalMu, 'actualMu', actualMu, ...
        'actualOuterFreeSurfaceRingRadius', actualRingRadius, ...
        'radii', radius, 'radialNormal', radialNormal, ...
        'tangentialNormal', tangentialNormal, 'verticalNormal', verticalNormal, ...
        'continuation', continuation, 'cacheSpecification', cacheSpecification);
end

function policy = outer_policy()
    policy = struct( ...
        'modelName', 'first-order-local-asymptotic-absorbing-bc-v1', ...
        'isExactDtN', false, ...
        'timeConvention', 'exp(+i*omega*t)', ...
        'outgoingRadialConvention', 'H_m^(2)(kappa*r)', ...
        'complexDispersionEquation', ...
        'kappa*tanh(kappa*h)=nu0*(1-i*mu_outer_actual)', ...
        'branchStart', 'positive-real-physical-k-at-mu=0', ...
        'branchConstraint', 'Re(kappa)>0;Im(kappa)<=1e-12', ...
        'failurePolicy', 'hard-stop-no-real-fallback-no-root-projection', ...
        'dispersionRelativeResidualTolerance', 1e-11, ...
        'newtonRelativeTolerance', 1e-13, ...
        'newtonMaximumIterationsPerStep', 40, ...
        'continuationBaseMuStep', 0.05, ...
        'continuationMaximumRefinements', 7, ...
        'muOuterDefinition', ...
        'nominal_mu*(outermost_FS_panel_center_ring_position)^2', ...
        'nominalMuDefinition', ...
        'resolved Phase2.2 actualMu0; legacy absent gives min(2.5,max(cfg.mu0,12/(k*width)))', ...
        'sourceLocalRadiusDefinition', 'r=hypot(source_center_x,source_center_y)', ...
        'storedNormalFormula', 'n_r=n dot e_r;n_theta=n dot e_theta;n_z=n dot e_z', ...
        'gammaFormula', ...
        'gamma=(-n_r)*(i*kappa+1/(2*r))+n_z*kappa*tanh(kappa*(z+h))', ...
        'canonicalSourceColumn', 'D+gamma_j*S', ...
        'radialDecayOrder', 1, ...
        'normalTangentialTolerance', 1e-10, ...
        'modeErrorReference', ...
        'exact-cylinder-gamma=-kappa*H_m^(2)''/H_m^(2),m=0,1');
end

function [kappa, info] = continue_complex_kappa(physicalK, depth, nu0, targetMu, policy)
    baseSteps = max(8, ceil(targetMu / policy.continuationBaseMuStep));
    success = false;
    for refinementLevel = 0:policy.continuationMaximumRefinements
        stepCount = baseSteps * 2 ^ refinementLevel;
        muPath = linspace(0, targetMu, stepCount + 1);
        candidate = physicalK;
        totalIterations = 0;
        maximumIterations = 0;
        minimumRealPart = real(candidate);
        maximumImagPart = imag(candidate);
        success = true;
        for stepIndex = 2:numel(muPath)
            rightHandSide = nu0 * (1 - 1i * muPath(stepIndex));
            [candidate, iterations, converged] = newton_complex_dispersion( ...
                candidate, depth, rightHandSide, policy);
            totalIterations = totalIterations + iterations;
            maximumIterations = max(maximumIterations, iterations);
            minimumRealPart = min(minimumRealPart, real(candidate));
            maximumImagPart = max(maximumImagPart, imag(candidate));
            relativeResidual = abs(candidate * tanh(candidate * depth) - ...
                rightHandSide) / max(abs(rightHandSide), eps);
            if ~converged || relativeResidual > ...
                    policy.dispersionRelativeResidualTolerance || ...
                    real(candidate) <= 0 || imag(candidate) > 1e-12
                success = false;
                break;
            end
        end
        if success, break; end
    end
    assert(success, 'CRESTU:OuterABCComplexDispersion', ...
        ['Complex-kappa continuation failed; no real-k fallback or root ', ...
        'projection is permitted.']);
    kappa = candidate;
    rightHandSide = nu0 * (1 - 1i * targetMu);
    relativeResidual = abs(kappa * tanh(kappa * depth) - rightHandSide) / ...
        max(abs(rightHandSide), eps);
    assert(real(kappa) > 0 && imag(kappa) <= 1e-12 && ...
        relativeResidual <= policy.dispersionRelativeResidualTolerance, ...
        'CRESTU:OuterABCBranch', 'Complex-kappa branch/residual contract failed.');
    info = struct('acceptedSteps', stepCount, ...
        'refinementLevel', refinementLevel, 'totalIterations', totalIterations, ...
        'maximumIterations', maximumIterations, ...
        'relativeResidual', relativeResidual, ...
        'minimumRealPart', minimumRealPart, 'maximumImagPart', maximumImagPart);
end

function [kappa, iterations, converged] = ...
        newton_complex_dispersion(initial, depth, rightHandSide, policy)
    kappa = initial;
    converged = false;
    for iterations = 1:policy.newtonMaximumIterationsPerStep
        kh = kappa * depth;
        functionValue = kappa * tanh(kh) - rightHandSide;
        derivative = tanh(kh) + kappa * depth / cosh(kh) ^ 2;
        if abs(derivative) <= 100 * eps, return; end
        update = functionValue / derivative;
        kappa = kappa - update;
        if abs(functionValue) / max(abs(rightHandSide), eps) <= ...
                policy.newtonRelativeTolerance && abs(update) <= ...
                policy.newtonRelativeTolerance * max(1, abs(kappa))
            converged = true;
            return;
        end
    end
end

function gamma = exact_cylinder_gamma(mode, kappa, radius)
    argument = kappa .* radius;
    hankel = besselh(mode, 2, argument);
    derivative = mode ./ argument .* hankel - besselh(mode + 1, 2, argument);
    gamma = -kappa .* derivative ./ hankel;
    assert(all(isfinite(real(gamma)) & isfinite(imag(gamma))), ...
        'CRESTU:OuterABCModeReference', 'Exact-cylinder reference is nonfinite.');
end

function hashText = numeric_sha256(values)
    values = reshape(double(values), 1, []);
    bytes = typecast(values, 'uint8');
    hashText = sha256_hash(bytes);
end
