function figure_handle = Plot_MeanDrift_Comparison(summary_file)
% PLOT_MEANDRIFT_COMPARISON Compare near/far drift and all Pinkster terms.
%
% Syntax:
%   figure_handle = Plot_MeanDrift_Comparison(summary_file)
%
% Inputs:
%   summary_file  : [char|string] Single-sphere convergence summary MAT-file.
%
% Outputs:
%   figure_handle : [Figure] Mean-drift validation and contribution figure.
%
% Mathematical Reference:
%   Pinkster near-field direct pressure and Maruo-Newman momentum conservation.
%% Stage 1: Initialize inputs and dependencies

    common_directory = fileparts(mfilename('fullpath'));

%% Stage 2: Run the core calculation

    case_wave_directory = fileparts(common_directory);
    code_root = fileparts(case_wave_directory);
    addpath(common_directory, fullfile(code_root,'5.Force'), fullfile(code_root,'6.MeanDriftLoads'));
    if nargin < 1 || isempty(summary_file)
        summary_file = fullfile(case_wave_directory,'Case1_SingleSphere_Convergence', ...
'SingleSphere_Convergence_Summary.mat');
    end
    loaded = load(summary_file,'summary');
    fine_result = loaded.summary.results{3};
    reference = loaded.summary.reference_irr3;
    omega = fine_result.omegas;
    heading_zero = find(abs(fine_result.headings) < 1.0e-12, 1);
    normalization = fine_result.drift.normalization;
    near_cd = squeeze(fine_result.drift.near(1, heading_zero, :)) / normalization;
    far_cd = squeeze(fine_result.drift.far(1, heading_zero, :)) / normalization;
    reference_index = match_reference_grid(reference.drift_pressure.omegas, omega);
    reference_heading = find(abs(reference.drift_pressure.headings) < 1.0e-12, 1);
    pressure_cd = squeeze(reference.drift_pressure.Cd(1, reference_heading, reference_index));
    momentum_cd = squeeze(reference.drift_momentum.Cd(1, reference_heading, reference_index));

    term_names = {'term_waterline','term_quadratic_velocity', ...
'term_rotation_force','term_translation_gradient'};
    term_values = zeros(numel(omega), 4);
    for frequency_index = 1:numel(omega)
        detail = fine_result.drift.near_details{frequency_index, 1};
        for term_index = 1:4
            term_values(frequency_index, term_index) = ...
                detail.(term_names{term_index})(1, heading_zero) / normalization;
        end
    end

    figure_handle = figure('Color','w','Position', [90, 90, 1350, 850], ...
'Name','Second-order mean-drift method comparison');
    layout = tiledlayout(2, 2,'TileSpacing','compact','Padding','compact');
    nexttile(layout);
    hold on;
    plot(omega, near_cd,'b-o','DisplayName','CRESTU near field');
    plot(omega, far_cd,'r-s','DisplayName','CRESTU far field');
    plot(omega, pressure_cd,'k--','DisplayName','WAMIT pressure method');
    plot(omega, momentum_cd,'m:','LineWidth', 1.7,'DisplayName','WAMIT momentum method');
    format_axis('Mean-drift method comparison','Mean surge drift coefficient C_d');

    nexttile(layout);
    plot(omega, term_values,'LineWidth', 1.2,'Marker','o');
    grid on;
    xlabel('\omega (rad/s)','FontWeight','bold');
    ylabel('Term C_d','FontWeight','bold');
    title('Four Pinkster near-field terms','FontWeight','bold');
    legend({'Relative waterline elevation','Velocity-squared term','Rotation and excitation coupling','Translation and pressure-gradient coupling'}, ...
'Location','best');

    nexttile(layout);
    near_error = 100 * abs(near_cd - pressure_cd) ./ max(abs(pressure_cd), 1.0e-10);
    far_error = 100 * abs(far_cd - momentum_cd) ./ max(abs(momentum_cd), 1.0e-10);
    semilogy(omega, near_error,'b-o', omega, far_error,'r-s');
    grid on;
    xlabel('\omega (rad/s)','FontWeight','bold');
    ylabel('Relative error (%)','FontWeight','bold');
    title('Relative error against WAMIT','FontWeight','bold');
    legend({'Near field / pressure method','Far field / momentum method'},'Location','best');

    nexttile(layout);
    consistency = 100 * abs(near_cd - far_cd) ./ max(0.5 * (abs(near_cd) + abs(far_cd)), 1.0e-10);
    plot(omega, consistency,'k-d','LineWidth', 1.2);
    grid on;
    xlabel('\omega (rad/s)','FontWeight','bold');
    ylabel('Near/far-field difference (%)','FontWeight','bold');
    title('Near/far-field consistency residual','FontWeight','bold');

    title(layout,'CRESTU second-order mean-drift validation','FontWeight','bold');
    output_file = fullfile(fileparts(summary_file),'Plot_MeanDrift_Comparison.png');
    export_chinese_figure(figure_handle, output_file);
    fprintf('[OK] Exported mean-drift comparison: %s\n', output_file);
    if nargout == 0
        clear figure_handle
    end
end

function format_axis(plot_title, y_label)
% FORMAT_AXIS Apply consistent Chinese mean-drift plot formatting.
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
