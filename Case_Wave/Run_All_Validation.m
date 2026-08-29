function validation = Run_All_Validation(force_recompute)
% RUN_ALL_VALIDATION Run the complete CRESTU validation suite.
%
% Syntax:
%   validation = Run_All_Validation(force_recompute)
%
% Description:
%   Runs single-sphere convergence and two-sphere interaction cases, then
%   exports CRESTU/WAMIT reconciliation tables.
%
% Inputs:
%   force_recompute - Recompute compatible potential solutions, logical scalar [-].
%
% Outputs:
%   validation      - Case results and validation metrics, structure [-].
%
% Governing Equations / Theory:
%   Uses linear radiation, diffraction, rigid-body response, and mean drift.
%
% References:
%   - WAMIT reference solutions supplied with the project.
%   - CRESTU theory and technical manual.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

    %% Stage 1: Apply defaults and validate suite directories

    if nargin < 1
        force_recompute = false;
    end

    validateattributes(force_recompute, {'logical', 'double'}, {'scalar'});
    caseWaveDirectory = fileparts(mfilename('fullpath'));
    singleCaseDirectory = fullfile(caseWaveDirectory, ...
        'Case1_SingleSphere_Convergence');
    multiCaseDirectory = fullfile(caseWaveDirectory, ...
        'Case2_TwoSpheres_Interaction');
    commonScriptDirectory = fullfile(caseWaveDirectory, 'Common_Scripts');
    requiredDirectories = {singleCaseDirectory, multiCaseDirectory, ...
        commonScriptDirectory};

    for directoryIndex = 1:numel(requiredDirectories)
        assert(isfolder(requiredDirectories{directoryIndex}), ...
            'CRESTU:MissingValidationDirectory', ...
            'Validation directory was not found: %s', ...
            requiredDirectories{directoryIndex});
    end

    addpath(singleCaseDirectory, multiCaseDirectory, commonScriptDirectory);

    %% Stage 2: Execute validation cases

    fprintf('[INFO] Start the complete CRESTU validation suite.\n');
    validation = struct();
    validation.single_body = Run_SingleSphere_Convergence(force_recompute);
    validation.two_body = Run_TwoSpheres_Interaction(force_recompute);

    %% Stage 3: Export and report reconciliation tables

    validation.metrics = Build_Validation_Tables();
    metricNames = ["A33"; "B33"; "Fz"; "HeaveRAO"];
    meanErrorPercent = validation.metrics.single_mean_error_pct(:); % [%]
    reconciliation = table(metricNames, meanErrorPercent, ...
        'VariableNames', {'Metric', 'MeanRelativeErrorPercent'});
    fprintf('[INFO] Validation reconciliation summary:\n');
    disp(reconciliation);
    fprintf('[OK] Complete validation suite finished.\n');
end
