function results = HydroMain(config_file)
% HYDROMAIN Run one CRESTU hydrodynamic case from a configuration file.
%
% Syntax:
%   results = HydroMain(config_file)
%
% Description:
%   Run one CRESTU hydrodynamic case from a configuration file.
%   The implementation preserves the CRESTU solver and data-field contracts.
%
% Inputs:
%   config_file - Path to a CRESTU configuration file, character vector or string scalar [-].
%
% Outputs:
%   results - Complete CRESTU frequency-domain result structure [-].
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
%HYDROMAIN Entry point for the CRESTU frequency-domain Rankine solver.
%% Stage 1: Initialize inputs and dependencies

    here = fileparts(mfilename('fullpath'));

%% Stage 2: Run the core calculation

    root = fileparts(here);
    addpath(fullfile(root,'1.Input'),fullfile(root,'2.Mesh'),fullfile(root,'3.HessSmith'), ...
        fullfile(root,'4.Potential'),fullfile(root,'5.Force'), ...
        fullfile(root,'6.MeanDriftLoads'),here);
    if nargin<1 || isempty(config_file)
        config_file = fullfile(here,'CRESTU.cfg');
    end
    results = run_frequency_domain_case(config_file);
end
