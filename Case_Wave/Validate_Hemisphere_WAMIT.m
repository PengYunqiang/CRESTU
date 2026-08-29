function validation = Validate_Hemisphere_WAMIT(run_solver,config_file)
% VALIDATE_HEMISPHERE_WAMIT Validate hemisphere hydrodynamic results against WAMIT data.
%
% Syntax:
%   validation = Validate_Hemisphere_WAMIT(run_solver,config_file)
%
% Description:
%   Validate hemisphere hydrodynamic results against WAMIT data.
%   The implementation preserves the CRESTU solver and data-field contracts.
%
% Inputs:
%   run_solver - Flag that enables a new solver run, logical scalar [-].
%   config_file - Path to a CRESTU configuration file, character vector or string scalar [-].
%
% Outputs:
%   validation - Validation results and error metrics structure [-].
%
% Governing Equations / Theory:
%   Uses the linear potential-flow, source-panel, or post-processing
%   formulation stated in the CRESTU theory and technical manual.
%
% References:
%   - CRESTU theory and technical manual.
%   - Standard boundary-element and marine hydrodynamics references.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)
%VALIDATE_HEMISPHERE_WAMIT Compare saved/new hemi_D10 results with WAMIT .1.
%% Stage 1: Initialize inputs and dependencies

    if nargin < 1
        run_solver = false;
    end

%% Stage 2: Run the core calculation
    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);
    addpath(fullfile(root,'1.Input'),fullfile(root,'2.Mesh'), ...
        fullfile(root,'4.Potential'),fullfile(root,'5.Force'),here);
    if nargin<2 || isempty(config_file)
        if run_solver
            config_file = fullfile(here,'CRESTU_Validation.cfg');
        else
            config_file = fullfile(here,'CRESTU.cfg');
        end
    end
    cfg = read_config(config_file);
    if run_solver
        results = HydroMain(cfg.config_file);
    elseif exist(cfg.files.results,'file')
        data = load(cfg.files.results,'results');
        results = data.results;
    else
        results = [];
    end
    wamit_file = fullfile(root,'..','WAMIT','FullSphereIRR0','hemisphere.1');
    reference = read_wamit_first_order(wamit_file,cfg.rho);
    validation = struct('reference',reference,'solver_results_available',~isempty(results), ...
'comparison',[],'sanity',struct());
    if isempty(results)
        warning('CRESTU:NoSolverResults','No hemi_D10 result cache exists; WAMIT reference ingestion only.');
        return;
    end
    validation.comparison = compare_wamit_results(results,reference);
    A33 = squeeze(results.added_mass(3,3,:));
    B33 = squeeze(results.damping(3,3,:));
    validation.sanity = struct('A33_positive',all(A33>0),'B33_positive',all(B33 >= 0), ...
'A_symmetry_max',max(arrayfun(@(k)norm(results.added_mass(:,:,k) - results.added_mass(:,:,k).','fro') / ...
            max(norm(results.added_mass(:,:,k),'fro'),eps),1:numel(results.omegas))), ...
'B_symmetry_max',max(arrayfun(@(k)norm(results.damping(:,:,k) - results.damping(:,:,k).','fro') / ...
            max(norm(results.damping(:,:,k),'fro'),eps),1:numel(results.omegas))), ...
'low_frequency_B33',B33(1),'high_frequency_B33',B33(end));
end
