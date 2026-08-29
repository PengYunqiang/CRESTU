function figure_handle = Plot_SingleBody_Results(summary_file)
% PLOT_SINGLEBODY_RESULTS Compare three CRESTU grids with WAMIT references.
%
% Syntax:
%   figure_handle = Plot_SingleBody_Results(summary_file)
%
% Inputs:
%   summary_file  : [char|string] SingleSphere_Convergence_Summary.mat path.
%
% Outputs:
%   figure_handle : [Figure] Six-panel comparison figure.
%
% Mathematical Reference:
%   Direct comparison of dimensional hydrodynamic coefficients and response amplitude operators.
%% Stage 1: Initialize inputs and dependencies

    [~, case_wave_directory, ~] = initialize_plot_paths();
    if nargin < 1 || isempty(summary_file)
        summary_file = fullfile(case_wave_directory,'Case1_SingleSphere_Convergence', ...
'SingleSphere_Convergence_Summary.mat');
    end

%% Stage 2: Run the core calculation
    loaded = load(summary_file,'summary');
    summary = loaded.summary;
    results = summary.results;
    omega = results{3}.omegas;
    headings = results{3}.headings;
    heading_zero = find(abs(headings) < 1.0e-12, 1);
    colors = lines(3);

    reference_zero = summary.reference_irr0;
    reference_three = summary.reference_irr3;
    reference_index_zero = match_reference_grid(reference_zero.radiation.omegas, omega);
    reference_index_three = match_reference_grid(reference_three.radiation.omegas, omega);
    reference_heading_zero = find(abs(reference_three.excitation.headings) < 1.0e-12, 1);

    figure_handle = figure('Color','w','Position', [80, 80, 1450, 820], ...
'Name','Single-sphere hydrodynamic comparison with WAMIT');
    layout = tiledlayout(2, 3,'TileSpacing','compact','Padding','compact');

    nexttile(layout);
    hold on;
    for level_index = 1:3
        plot(omega, squeeze(results{level_index}.added_mass(3, 3, :)),'-o', ...
'Color', colors(level_index, :),'DisplayName', ['CRESTU-' summary.levels{level_index}]);
    end
    plot(omega, squeeze(reference_zero.radiation.added_mass(3, 3, reference_index_zero)),'k--', ...
'LineWidth', 1.5,'DisplayName','WAMIT-IRR0');
    plot(omega, squeeze(reference_three.radiation.added_mass(3, 3, reference_index_three)),'m:', ...
'LineWidth', 1.8,'DisplayName','WAMIT-IRR3');
    format_axis('Heave added mass A_{33}','A_{33} (kg)');

    nexttile(layout);
    hold on;
    for level_index = 1:3
        plot(omega, squeeze(results{level_index}.damping(3, 3, :)),'-o', ...
'Color', colors(level_index, :),'DisplayName', ['CRESTU-' summary.levels{level_index}]);
    end
    plot(omega, squeeze(reference_zero.radiation.damping(3, 3, reference_index_zero)),'k--', ...
'LineWidth', 1.5,'DisplayName','WAMIT-IRR0');
    plot(omega, squeeze(reference_three.radiation.damping(3, 3, reference_index_three)),'m:', ...
'LineWidth', 1.8,'DisplayName','WAMIT-IRR3');
    format_axis('Heave radiation damping B_{33}','B_{33} (kg/s)');

    nexttile(layout);
    hold on;
    fine_excitation = squeeze(abs(results{3}.excitation(3, heading_zero, :)));
    wamit_excitation = squeeze(abs(reference_three.excitation.force( ...
        3, reference_heading_zero, reference_index_three)));
    plot(omega, fine_excitation,'b-o','DisplayName','CRESTU-Fine');
    plot(omega, wamit_excitation,'k--','LineWidth', 1.5,'DisplayName','WAMIT-IRR3');
    format_axis('Heave excitation at zero-degree heading','|F_{z,exc}| (N)');

    nexttile(layout);
    hold on;
    fine_rao = squeeze(results{3}.rao.amplitude(3, heading_zero, :));
    wamit_rao = squeeze(reference_three.rao.amplitude(3, reference_heading_zero, reference_index_three));
    plot(omega, fine_rao,'b-o','DisplayName','CRESTU-Fine');
    plot(omega, wamit_rao,'k--','LineWidth', 1.5,'DisplayName','WAMIT-IRR3');
    format_axis('Heave response amplitude operator','|\Xi_3| (m/m)');

    nexttile(layout);
    hold on;
    normalization = results{3}.drift.normalization;
    near_cd = squeeze(results{3}.drift.near(1, heading_zero, :)) / normalization;
    far_cd = squeeze(results{3}.drift.far(1, heading_zero, :)) / normalization;
    drift_index = match_reference_grid(reference_three.drift_pressure.omegas, omega);
    drift_heading_zero = find(abs(reference_three.drift_pressure.headings) < 1.0e-12, 1);
    wamit_near = squeeze(reference_three.drift_pressure.Cd(1, drift_heading_zero, drift_index));
    wamit_far = squeeze(reference_three.drift_momentum.Cd(1, drift_heading_zero, drift_index));
    plot(omega, near_cd,'b-o','DisplayName','CRESTU near field');
    plot(omega, far_cd,'r-s','DisplayName','CRESTU far field');
    plot(omega, wamit_near,'k--','DisplayName','WAMIT pressure method');
    plot(omega, wamit_far,'m:','LineWidth', 1.5,'DisplayName','WAMIT momentum method');
    format_axis('Mean surge drift coefficient','C_d (-)');

    nexttile(layout);
    hold on;
    fine_a = squeeze(results{3}.added_mass(3, 3, :));
    fine_b = squeeze(results{3}.damping(3, 3, :));
    reference_a = squeeze(reference_three.radiation.added_mass(3, 3, reference_index_three));
    reference_b = squeeze(reference_three.radiation.damping(3, 3, reference_index_three));
    error_a = 100 * abs(fine_a - reference_a) ./ max(abs(reference_a), eps);
    error_b = 100 * abs(fine_b - reference_b) ./ max(abs(reference_b), eps);
    error_f = 100 * abs(fine_excitation - wamit_excitation) ./ max(abs(wamit_excitation), eps);
    error_r = 100 * abs(fine_rao - wamit_rao) ./ max(abs(wamit_rao), eps);
    plot(omega, error_a,'-o','DisplayName','A_{33}');
    plot(omega, error_b,'-s','DisplayName','B_{33}');
    plot(omega, error_f,'-^','DisplayName','|F_{z,exc}|');
    plot(omega, error_r,'-d','DisplayName','Heave RAO');
    format_axis('Fine-gridRelative error','Relative error (%)');

    title(layout,'CRESTU single-sphere grid convergence and WAMIT comparison (\omega=0.5:0.1:2.0 rad/s)', ...
'FontWeight','bold');
    output_file = fullfile(fileparts(summary_file),'Plot_SingleBody_Results.png');
    export_chinese_figure(figure_handle, output_file);
    fprintf('[OK] Exported single-sphere comparison: %s\n', output_file);

    if nargout == 0
        clear figure_handle
    end
end

function [common_directory, case_wave_directory, code_root] = initialize_plot_paths()
% INITIALIZE_PLOT_PATHS Resolve and add the standalone plotting dependencies.
%
% Syntax:
%   [common_directory, case_wave_directory, code_root] = initialize_plot_paths()
%
% Inputs:
%   None.
%
% Outputs:
%   common_directory    : [char] Common plotting-script directory.
%   case_wave_directory : [char] Case_Wave directory.
%   code_root           : [char] CRESTU Code directory.
%
% Mathematical Reference:
%   Path utility; no mathematical model is used.
%% Stage 1: Initialize inputs and dependencies

    common_directory = fileparts(mfilename('fullpath'));

%% Stage 2: Run the core calculation

    case_wave_directory = fileparts(common_directory);
    code_root = fileparts(case_wave_directory);
    addpath(common_directory, fullfile(code_root,'5.Force'), fullfile(code_root,'6.MeanDriftLoads'));
end

function format_axis(plot_title, y_label)
% FORMAT_AXIS Apply consistent Chinese plot formatting.
%
% Syntax:
%   format_axis(plot_title, y_label)
%
% Inputs:
%   plot_title : [char|string] Chinese subplot title.
%   y_label    : [char|string] Chinese vertical-axis label.
%
% Outputs:
%   None; the current axes are formatted.
%
% Mathematical Reference:
%   Visualization utility; no mathematical model is used.
%% Stage 1: Initialize inputs and dependencies

    grid on;
    box on;
    xlabel('\omega (rad/s)','FontWeight','bold');
    ylabel(y_label,'FontWeight','bold');
    title(plot_title,'FontWeight','bold');
    legend('Location','best');
end
