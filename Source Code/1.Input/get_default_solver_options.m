function options = get_default_solver_options()
% GET_DEFAULT_SOLVER_OPTIONS Central defaults for production and validation controls.
% Defaults reproduce the pre-parameterization implementation.  Theory
% constants (pi, Green factors, and matrix dimensions) remain in algorithms.

    options = struct();

    % Controls that can change generated geometry or production boundary data.
    options.numerics = struct();
    options.numerics.mesh = struct();
    options.numerics.mesh.mode = 'explicit';
    options.numerics.mesh.targetPanelsPerWavelength = 24;
    options.numerics.mesh.minFreeSurfaceRadialLayers = 3;
    options.numerics.mesh.maxFreeSurfaceRadialLayers = 256;
    options.numerics.mesh.minBottomRadialLayers = 3;
    options.numerics.mesh.maxBottomRadialLayers = 256;
    options.numerics.mesh.minThetaPanels = 4;
    options.numerics.mesh.maxThetaPanels = 4096;
    options.numerics.mesh.freeSurfaceNearGrowthRatio = 1.0;
    options.numerics.mesh.freeSurfaceSmoothingMaxIterations = 250;
    options.numerics.mesh.freeSurfaceSmoothingRelaxation = 1.15;
    options.numerics.mesh.freeSurfaceSmoothingDenominatorTolerance = 1e-10;
    options.numerics.mesh.seabed = struct();
    options.numerics.mesh.seabed.waterlineLengthTolerance = 1e-12;
    options.numerics.mesh.seabed.cornerAngleThresholdDeg = 30;
    options.numerics.mesh.seabed.centerBlend = 0.5;
    options.numerics.mesh.seabed.radialLayerFraction = 0.25;
    options.numerics.mesh.seabed.minimumRadialLayers = 3;
    options.numerics.mesh.seabed.smoothingIterations = 25;
    % 1.0 reproduces the original full Jacobi replacement; lower values
    % provide an explicitly configurable under-relaxation.
    options.numerics.mesh.seabed.smoothingWeight = 1.0;
    options.numerics.mesh.seabed.segmentLengthTolerance = 1e-8;
    options.numerics.mesh.reducedSeabedRadialExponent = 0.8;
    options.numerics.mesh.farfieldDefaultVerticalLayers = 8;
    options.numerics.mesh.infiniteDepthBoundaryDepthFactor = 0.8;
    options.numerics.mesh.resolutionToleranceFraction = 0.05;
    options.numerics.mesh.outerMaximumCandidates = 256;
    options.numerics.mesh.outerRingToleranceScale = 2e-8;
    options.numerics.mesh.outerPanelDegeneracyTolerance = 1e-12;
    options.numerics.mesh.outerPanelPlanarityTolerance = 1e-13;
    options.numerics.mesh.meshMergeTolerance = 1e-10;
    options.numerics.mesh.panelGeometryTolerance = 1e-12;
    options.numerics.mesh.symmetryTolerance = 1e-9;
    options.numerics.mesh.diskOuterRingToleranceScale = 1e-6;
    options.numerics.mesh.gridRadialDivisor = 14;
    options.numerics.mesh.gridMinimumSpacingM = 0.75;
    options.numerics.mesh.interiorDomainFraction = 0.985;
    options.numerics.mesh.multibodyMinimumOuterThetaCount = 48;
    options.numerics.mesh.diskFallbackThetaCount = 64;
    options.numerics.mesh.waterlineHorizontalNormalThreshold = 0.2;
    options.numerics.mesh.waterlineNodeDistanceTolerance = 1e-3;
    options.numerics.mesh.outerTruncationTopRadiusPerWavelength = 1.5;
    options.numerics.mesh.outerTruncationBottomRadiusPerWavelength = 1.5;
    options.numerics.mesh.outerTruncationMinimumRadiusPerWavelength = 1.5;
    options.numerics.mesh.spongeMinimumWidthWavelengths = 1.5;
    options.numerics.mesh.spongeAutoMuMaximum = 2.5;
    options.numerics.mesh.spongeAutoMuTargetIntegral = 12;
    options.numerics.mesh.spongeTargetAttenuationExponent = 4;
    options.numerics.mesh.spongeRadiusExpansionRelativeTolerance = 1e-12;

    options.numerics.assembly = struct( ...
        'rankineFieldBlockSize', 1024, ...
        'useFusedSymmetry', false);

    options.numerics.solver = struct( ...
        'dispersionMaximumIterations', 50, ...
        'dispersionRelativeTolerance', 1e-12, ...
        'raoRcondWarningThreshold', 1e-12);

    options.numerics.outerABC = struct();
    options.numerics.outerABC.dispersionRelativeResidualTolerance = 1e-11;
    options.numerics.outerABC.newtonRelativeTolerance = 1e-13;
    options.numerics.outerABC.newtonMaximumIterationsPerStep = 40;
    options.numerics.outerABC.continuationBaseMuStep = 0.05;
    options.numerics.outerABC.continuationMaximumRefinements = 7;
    options.numerics.outerABC.normalTangentialTolerance = 1e-10;
    options.numerics.outerABC.derivativeBreakdownToleranceMultiplier = 100;
    options.numerics.outerABC.branchImaginaryTolerance = 1e-12;

    % Validation/diagnostic controls are deliberately separate from physics.
    options.validation = struct();
    options.validation.haskindThetaCount = 72;
    options.validation.physicalBoundaryEpsilonFactor = 0.35;
    options.validation.physicalBoundaryResidualThreshold = 0.25;
    options.validation.haskindPSDRelativeTolerance = 1e-10;
    options.validation.frequencyMatchToleranceRadPerS = 2e-4;
    options.validation.radiationFlux = struct( ...
        'radiusFractions', [0.3, 0.6, 0.9], ...
        'thetaLevels', [16, 24, 32], 'zLevels', [8, 10, 12], ...
        'radialStepWavelengthFraction', 0.02, ...
        'radialStepGapFraction', 0.2, ...
        'validationTolerance', 0.05, 'hermitianTolerance', 0.05);
    options.validation.resolutionToleranceFraction = 0.05;

    % Cache matching policy.  Schema/version values are not user physics.
    options.cache = struct('frequencyToleranceMultiplier', 8, ...
        'geometryToleranceMultiplier', 64);
end
