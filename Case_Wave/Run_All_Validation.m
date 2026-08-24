function validation = Run_All_Validation(force_recompute)
% RUN_ALL_VALIDATION Run the complete single- and two-body validation suite.
%
% Syntax:
%   validation = Run_All_Validation(force_recompute)
%
% Inputs:
%   force_recompute : [logical scalar] Recompute potentials instead of loading compatible caches.
%
% Outputs:
%   validation       : [struct] Single-body, two-body, and tabulated validation results.
%
% Mathematical Reference:
%   Frequency-domain radiation, diffraction, response, and mean-drift benchmarks against WAMIT.

% ========================================== %
% Initialize validation paths
% ========================================== %
    if nargin < 1
        force_recompute = false;
    end
    case_wave_directory = fileparts(mfilename('fullpath'));
    single_directory = fullfile(case_wave_directory, 'Case1_SingleSphere_Convergence');
    multi_directory = fullfile(case_wave_directory, 'Case2_TwoSpheres_Interaction');
    common_directory = fullfile(case_wave_directory, 'Common_Scripts');
    addpath(single_directory, multi_directory, common_directory);

% ========================================== %
% Execute cases and export comparison tables
% ========================================== %
    validation = struct();
    validation.single_body = Run_SingleSphere_Convergence(force_recompute);
    validation.two_body = Run_TwoSpheres_Interaction(force_recompute);
    validation.metrics = Build_Validation_Tables();
end
