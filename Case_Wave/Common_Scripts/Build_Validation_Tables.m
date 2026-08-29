function metrics = Build_Validation_Tables(single_summary_file, multi_summary_file)
% BUILD_VALIDATION_TABLES Export per-frequency CRESTU/WAMIT error tables.
%
% Syntax:
%   metrics = Build_Validation_Tables(single_summary_file, multi_summary_file)
%
% Inputs:
%   single_summary_file : [char|string] Single-sphere convergence summary MAT-file.
%   multi_summary_file  : [char|string] Two-sphere interaction summary MAT-file.
%
% Outputs:
%   metrics             : [struct] Aggregate errors, symmetry residuals, and CSV paths.
%
% Mathematical Reference:
%   Relative error = 100*abs(CRESTU-WAMIT)/max(abs(WAMIT),eps).
%% Stage 1: Resolve paths and apply defaults

    common_directory = fileparts(mfilename('fullpath'));

%% Stage 2: Run the core calculation

    case_wave_directory = fileparts(common_directory);
    code_root = fileparts(case_wave_directory);
    addpath(common_directory, fullfile(code_root,'5.Force'), fullfile(code_root,'6.MeanDriftLoads'));
    if nargin < 1 || isempty(single_summary_file)
        single_summary_file = fullfile(case_wave_directory,'Case1_SingleSphere_Convergence', ...
'SingleSphere_Convergence_Summary.mat');
    end
    if nargin < 2 || isempty(multi_summary_file)
        multi_summary_file = fullfile(case_wave_directory,'Case2_TwoSpheres_Interaction', ...
'TwoSpheres_Interaction_Summary.mat');
    end
    assert(isfile(single_summary_file), 'CRESTU:MissingValidationSummary', ...
        'Single-sphere summary was not found: %s', single_summary_file);
    assert(isfile(multi_summary_file), 'CRESTU:MissingValidationSummary', ...
        'Two-sphere summary was not found: %s', multi_summary_file);
    fprintf('[INFO] Build validation tables from completed case summaries.\n');

%% Stage 2: Build the single-sphere frequency table
    single_data = load(single_summary_file,'summary');
    single = single_data.summary;
    fine = single.results{3};
    medium = single.results{2};
    reference = single.reference_irr3;
    frequency_index = match_reference_grid(reference.radiation.omegas, fine.omegas);
    drift_index = match_reference_grid(reference.drift_pressure.omegas, fine.omegas);
    omega = fine.omegas(:);

    added_mass = squeeze(fine.added_mass(3, 3, :));
    added_mass_reference = squeeze(reference.radiation.added_mass(3, 3, frequency_index));
    damping = squeeze(fine.damping(3, 3, :));
    damping_reference = squeeze(reference.radiation.damping(3, 3, frequency_index));
    excitation = squeeze(abs(fine.excitation(3, 1, :)));
    excitation_reference = squeeze(abs(reference.excitation.force(3, 1, frequency_index)));
    heave_rao = squeeze(fine.rao.amplitude(3, 1, :));
    heave_rao_reference = squeeze(reference.rao.amplitude(3, 1, frequency_index));
    drift_near_cd = squeeze(fine.drift.near(1, 1, :)) / fine.drift.normalization;
    drift_far_cd = squeeze(fine.drift.far(1, 1, :)) / fine.drift.normalization;
    drift_pressure_cd = squeeze(reference.drift_pressure.Cd(1, 1, drift_index));
    drift_momentum_cd = squeeze(reference.drift_momentum.Cd(1, 1, drift_index));

    single_table = table(omega, added_mass, added_mass_reference, ...
        relative_error(added_mass, added_mass_reference), damping, damping_reference, ...
        relative_error(damping, damping_reference), excitation, excitation_reference, ...
        relative_error(excitation, excitation_reference), heave_rao, heave_rao_reference, ...
        relative_error(heave_rao, heave_rao_reference), drift_near_cd, drift_pressure_cd, ...
        drift_far_cd, drift_momentum_cd, ...
        relative_error(squeeze(fine.added_mass(3, 3, :)), squeeze(medium.added_mass(3, 3, :))), ...
        relative_error(squeeze(fine.damping(3, 3, :)), squeeze(medium.damping(3, 3, :))), ...
'VariableNames', {'omega_rad_s','A33_CRESTU','A33_WAMIT','A33_error_pct', ...
'B33_CRESTU','B33_WAMIT','B33_error_pct','Fz_CRESTU_N','Fz_WAMIT_N', ...
'Fz_error_pct','Heave_RAO_CRESTU','Heave_RAO_WAMIT','Heave_RAO_error_pct', ...
'Drift_near_Cd','Drift_WAMIT_pressure_Cd','Drift_far_Cd', ...
'Drift_WAMIT_momentum_Cd','A33_fine_medium_pct','B33_fine_medium_pct'});
    single_csv = fullfile(fileparts(single_summary_file),'SingleSphere_Validation_Metrics.csv');
    writetable(single_table, single_csv);

%% Stage 3: Build the two-sphere coefficient and RAO table
    multi_data = load(multi_summary_file,'summary');
    multi = multi_data.summary;
    multi_result = multi.results;
    multi_reference = multi.reference;
    multi_index = match_reference_grid(multi_reference.radiation.omegas, multi_result.omegas);
    omega_multi = multi_result.omegas(:);
    coefficient_pairs = [1, 1; 3, 3; 1, 7; 3, 9];
    coefficient_names = {'11','33','17','39'};
    multi_table = table(omega_multi,'VariableNames', {'omega_rad_s'});
    for pair_index = 1:size(coefficient_pairs, 1)
        row = coefficient_pairs(pair_index, 1);
        column = coefficient_pairs(pair_index, 2);
        suffix = coefficient_names{pair_index};
        value_a = squeeze(multi_result.added_mass(row, column, :));
        reference_a = squeeze(multi_reference.radiation.added_mass(row, column, multi_index));
        value_b = squeeze(multi_result.damping(row, column, :));
        reference_b = squeeze(multi_reference.radiation.damping(row, column, multi_index));
        multi_table.(['A' suffix '_CRESTU']) = value_a;
        multi_table.(['A' suffix '_WAMIT']) = reference_a;
        multi_table.(['A' suffix '_error_pct']) = relative_error(value_a, reference_a);
        multi_table.(['B' suffix '_CRESTU']) = value_b;
        multi_table.(['B' suffix '_WAMIT']) = reference_b;
        multi_table.(['B' suffix '_error_pct']) = relative_error(value_b, reference_b);
    end
    heading_pairs = [0, 0; 90, 90];
    for heading_index = 1:size(heading_pairs, 1)
        crestu_heading = find(abs(multi_result.headings - heading_pairs(heading_index, 1)) < 1.0e-12, 1);
        wamit_heading = find(abs(multi_reference.rao.headings - heading_pairs(heading_index, 2)) < 1.0e-12, 1);
        crestu_rao = squeeze(multi_result.rao.amplitude(:, crestu_heading, :));
        wamit_rao = squeeze(multi_reference.rao.amplitude(:, wamit_heading, multi_index));
        nrmse = zeros(numel(omega_multi), 1);
        for frequency = 1:numel(omega_multi)
            nrmse(frequency) = 100 * norm(crestu_rao(:, frequency) - wamit_rao(:, frequency)) ...
 / max(norm(wamit_rao(:, frequency)), eps);
        end
        multi_table.(sprintf('RAO_NRMSE_%ddeg_pct', heading_pairs(heading_index, 1))) = nrmse;
    end
    multi_csv = fullfile(fileparts(multi_summary_file),'TwoSpheres_Validation_Metrics.csv');
    writetable(multi_table, multi_csv);

%% Stage 4: Aggregate metrics and export reconciliation files
    metrics = struct();
    metrics.single_csv = single_csv;
    metrics.multi_csv = multi_csv;
    metrics.single_mean_error_pct = [mean(single_table.A33_error_pct), ...
        mean(single_table.B33_error_pct), mean(single_table.Fz_error_pct), ...
        mean(single_table.Heave_RAO_error_pct)];
    metrics.single_omega1_error_pct = [single_table.A33_error_pct(6), ...
        single_table.B33_error_pct(6), single_table.Fz_error_pct(6), ...
        single_table.Heave_RAO_error_pct(6)];
    metrics.multi_median_error_pct = zeros(4, 2);
    for pair_index = 1:4
        suffix = coefficient_names{pair_index};
        metrics.multi_median_error_pct(pair_index, 1) = median(multi_table.(['A' suffix '_error_pct']));
        metrics.multi_median_error_pct(pair_index, 2) = median(multi_table.(['B' suffix '_error_pct']));
    end
    metrics.max_added_mass_symmetry_residual = max(arrayfun(@(index) ...
        norm(multi_result.added_mass(:, :, index) - multi_result.added_mass(:, :, index).','fro') ...
 / max(norm(multi_result.added_mass(:, :, index),'fro'), eps), 1:numel(omega_multi)));
    metrics.max_damping_symmetry_residual = max(arrayfun(@(index) ...
        norm(multi_result.damping(:, :, index) - multi_result.damping(:, :, index).','fro') ...
 / max(norm(multi_result.damping(:, :, index),'fro'), eps), 1:numel(omega_multi)));
    metrics_file = fullfile(case_wave_directory,'Validation_Metrics.mat');
    save(metrics_file,'metrics','single_table','multi_table');
    reconciliation = table(["A33"; "B33"; "Fz"; "HeaveRAO"], ...
        metrics.single_mean_error_pct(:), ...
        'VariableNames', {'Metric', 'MeanRelativeErrorPercent'});
    disp(reconciliation);
    fprintf('[OK] Validation tables exported:\n    %s\n    %s\n', single_csv, multi_csv);
end

function error_pct = relative_error(value, reference)
% RELATIVE_ERROR Compute elementwise absolute percentage error.
%
% Syntax:
%   error_pct = relative_error(value, reference)
%
% Inputs:
%   value     : [N x 1] CRESTU values.
%   reference : [N x 1] WAMIT values.
%
% Outputs:
%   error_pct : [N x 1] Absolute relative errors, in percent.
%
% Mathematical Reference:
%   error_pct = 100*abs(value-reference)/max(abs(reference),eps).
%% Stage 1: Initialize inputs and dependencies

    error_pct = 100 * abs(value - reference) ./ max(abs(reference), eps);
end
