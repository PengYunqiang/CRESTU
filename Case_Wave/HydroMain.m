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

    caseWaveDirectory = fileparts(mfilename('fullpath'));

%% 阶段 2: 加载真实源码目录并运行主求解

    projectDirectory = fileparts(caseWaveDirectory);
    sourceDirectory = fullfile(projectDirectory, 'Source Code');
    assert(isfolder(sourceDirectory), 'CRESTU:SourceDirectoryMissing', ...
        'CRESTU source directory was not found: %s', sourceDirectory);
    addpath(fullfile(sourceDirectory,'0.Tools Code'), ...
        fullfile(sourceDirectory,'1.Input'), ...
        fullfile(sourceDirectory,'2.Mesh'), ...
        fullfile(sourceDirectory,'3.HessSmith'), ...
        fullfile(sourceDirectory,'4.Potential'), ...
        fullfile(sourceDirectory,'5.Force'), ...
        fullfile(sourceDirectory,'6.MeanDriftLoads'),caseWaveDirectory);
    if nargin<1 || isempty(config_file)
        config_file = fullfile(caseWaveDirectory,'CRESTU.cfg');
    end
    results = run_frequency_domain_case(config_file);
end
