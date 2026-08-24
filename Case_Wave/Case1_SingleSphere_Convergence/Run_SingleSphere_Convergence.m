function summary = Run_SingleSphere_Convergence(force_recompute)
% RUN_SINGLESPHERE_CONVERGENCE Run the 16-frequency three-grid hemisphere study.
%
% Syntax:
%   summary = Run_SingleSphere_Convergence(force_recompute)
%
% Inputs:
%   force_recompute : [logical scalar] Ignore compatible result caches when true.
%
% Outputs:
%   summary         : [struct] Three CRESTU results, WAMIT references, runtimes, and mesh statistics.
%
% Mathematical Reference:
%   Grid-convergence assessment against WAMIT FullSphereIRR0 and FullSphereIRR3.
    if nargin < 1
        force_recompute = false;
    end
    [case_directory, ~, wamit_root] = initialize_paths();
    level_names = {'Coarse', 'Medium', 'Fine'};
    subdirectories = {'Mesh_Coarse', 'Mesh_Medium', 'Mesh_Fine'};
    config_names = {'Case1_Coarse.cfg', 'Case1_Medium.cfg', 'Case1_Fine.cfg'};
    result_set = cell(3, 1);
    runtime_seconds = zeros(3, 1);

    for level_index = 1:3
        config_file = fullfile(case_directory, subdirectories{level_index}, config_names{level_index});
        cfg = read_config(config_file);
        timer_id = tic;
        if ~force_recompute && is_compatible_result(cfg.files.results, cfg)
            loaded = load(cfg.files.results, 'results');
            result_set{level_index} = loaded.results;
        else
            result_set{level_index} = run_frequency_domain_case(config_file);
        end
        runtime_seconds(level_index) = toc(timer_id);
    end

    reference_irr0 = read_wamit_dataset(fullfile(wamit_root, 'FullSphereIRR0'), 10);
    reference_irr3 = read_wamit_dataset(fullfile(wamit_root, 'FullSphereIRR3'), 10);
    panel_counts = cellfun(@(item) item.stats.total_dofs, result_set);
    summary = struct('schema_version', 1, 'levels', {level_names}, 'results', {result_set}, ...
        'runtime_seconds', runtime_seconds, 'panel_counts', panel_counts, ...
        'reference_irr0', reference_irr0, 'reference_irr3', reference_irr3, ...
        'frequency_grid', 0.5:0.1:2.0, 'headings', [0, 45, 90]);
    summary_file = fullfile(case_directory, 'SingleSphere_Convergence_Summary.mat');
    save(summary_file, 'summary', '-v7.3');
    Plot_SingleBody_Results(summary_file);
    Plot_MeanDrift_Comparison(summary_file);
    fprintf('>>> Single-sphere convergence suite complete: %s\n', summary_file);
end

function [case_directory, code_root, wamit_root] = initialize_paths()
% INITIALIZE_PATHS Add CRESTU modules and resolve the read-only benchmark root.
%
% Syntax:
%   [case_directory, code_root, wamit_root] = initialize_paths()
%
% Inputs:
%   None.
%
% Outputs:
%   case_directory : [char] Directory containing this convergence suite.
%   code_root      : [char] CRESTU Code directory.
%   wamit_root     : [char] Read-only sibling WAMIT directory.
%
% Mathematical Reference:
%   Path utility; no mathematical model is used.
    case_directory = fileparts(mfilename('fullpath'));
    case_wave_directory = fileparts(case_directory);
    code_root = fileparts(case_wave_directory);
    wamit_root = fullfile(fileparts(code_root), 'WAMIT');
    addpath(fullfile(code_root, '1.Input'), fullfile(code_root, '2.Mesh'), ...
        fullfile(code_root, '3.HessSmith'), fullfile(code_root, '4.Potential'), ...
        fullfile(code_root, '5.Force'), fullfile(code_root, '6.MeanDriftLoads'), ...
        case_wave_directory, fullfile(case_wave_directory, 'Common_Scripts'));
end

function compatible = is_compatible_result(filename, cfg)
% IS_COMPATIBLE_RESULT Check whether a saved result matches the requested grid.
%
% Syntax:
%   compatible = is_compatible_result(filename, cfg)
%
% Inputs:
%   filename : [char] Result MAT-file path.
%   cfg      : [struct] Requested configuration.
%
% Outputs:
%   compatible : [logical scalar] True for a reusable result.
%
% Mathematical Reference:
%   Exact metadata validation with floating-point tolerance.
    compatible = false;
    if ~isfile(filename)
        return
    end
    loaded = load(filename, 'results');
    if ~isfield(loaded, 'results') || ~isfield(loaded.results, 'schema_version')
        return
    end
    compatible = loaded.results.schema_version >= 5 ...
        && isequal(size(loaded.results.omegas), size(cfg.freq.omegas)) ...
        && max(abs(loaded.results.omegas - cfg.freq.omegas)) < 1.0e-12 ...
        && isequal(loaded.results.headings, cfg.wave.headings) ...
        && isfield(loaded.results, 'drift') && loaded.results.drift.enabled;
end
