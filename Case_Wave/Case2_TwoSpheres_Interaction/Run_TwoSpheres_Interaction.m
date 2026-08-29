function summary = Run_TwoSpheres_Interaction(force_recompute)
% RUN_TWOSPHERES_INTERACTION Run and benchmark the 12-DOF two-sphere case.
%
% Syntax:
%   summary = Run_TwoSpheres_Interaction(force_recompute)
%
% Inputs:
%   force_recompute : [logical scalar] Ignore a compatible local result when true.
%
% Outputs:
%   summary         : [struct] CRESTU result, WAMIT IRR3 reference, runtime, and dimensions.
%
% Mathematical Reference:
%   Linear multi-body radiation/diffraction interaction with 6 DOF per body.
%% Stage 1: Initialize inputs and dependencies

    if nargin < 1
        force_recompute = false;
    end

%% Stage 2: Run the core calculation
    case_directory = fileparts(mfilename('fullpath'));
    case_wave_directory = fileparts(case_directory);
    code_root = fileparts(case_wave_directory);
    wamit_root = fullfile(fileparts(code_root),'WAMIT');
    addpath(fullfile(code_root,'1.Input'), fullfile(code_root,'2.Mesh'), ...
        fullfile(code_root,'3.HessSmith'), fullfile(code_root,'4.Potential'), ...
        fullfile(code_root,'5.Force'), fullfile(code_root,'6.MeanDriftLoads'), ...
        case_wave_directory, fullfile(case_wave_directory,'Common_Scripts'));

    config_file = fullfile(case_directory,'Case2_TwoSpheres.cfg');
    cfg = read_config(config_file);
    timer_id = tic;
    if ~force_recompute && isfile(cfg.files.results)
        loaded = load(cfg.files.results,'results');
        compatible = isfield(loaded,'results') && loaded.results.schema_version >= 5 ...
 && isequal(loaded.results.omegas, cfg.freq.omegas) ...
 && isequal(loaded.results.headings, cfg.wave.headings) ...
 && size(loaded.results.added_mass, 1) == 12 ...
 && loaded.results.drift.enabled;
    else
        compatible = false;
    end
    if compatible
        results = loaded.results;
    else
        results = run_frequency_domain_case(config_file);
    end
    runtime_seconds = toc(timer_id);

    reference = read_wamit_dataset(fullfile(wamit_root,'TwoSphereIRR3'), 10);
    summary = struct('schema_version', 1,'results', results,'reference', reference, ...
'runtime_seconds', runtime_seconds,'frequency_grid', 0.5:0.1:2.0, ...
'headings', [0, 45, 90],'dof_count', size(results.added_mass, 1));
    summary_file = fullfile(case_directory,'TwoSpheres_Interaction_Summary.mat');
    save(summary_file,'summary','-v7.3');
    Plot_MultiBody_Results(summary_file);
    fprintf('[OK] Two-sphere interaction suite completed: %s\n', summary_file);
end
