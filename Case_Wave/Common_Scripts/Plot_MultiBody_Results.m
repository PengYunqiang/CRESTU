function figure_handle = Plot_MultiBody_Results(summary_file)
% PLOT_MULTIBODY_RESULTS Compare 12-DOF two-sphere CRESTU and WAMIT results.
%
% Syntax:
%   figure_handle = Plot_MultiBody_Results(summary_file)
%
% Inputs:
%   summary_file  : [char|string] TwoSpheres_Interaction_Summary.mat path.
%
% Outputs:
%   figure_handle : [Figure] Self/cross coefficient and RAO comparison figure.
%
% Mathematical Reference:
%   Multi-body radiation interaction coefficients and coupled motion RAOs.
    common_directory = fileparts(mfilename('fullpath'));
    case_wave_directory = fileparts(common_directory);
    code_root = fileparts(case_wave_directory);
    addpath(common_directory, fullfile(code_root, '5.Force'), fullfile(code_root, '6.MeanDriftLoads'));
    if nargin < 1 || isempty(summary_file)
        summary_file = fullfile(case_wave_directory, 'Case2_TwoSpheres_Interaction', ...
            'TwoSpheres_Interaction_Summary.mat');
    end
    loaded = load(summary_file, 'summary');
    results = loaded.summary.results;
    reference = loaded.summary.reference;
    omega = results.omegas;
    reference_index = match_reference_grid(reference.radiation.omegas, omega);
    heading_zero = find(abs(results.headings) < 1.0e-12, 1);
    reference_heading_zero = find(abs(reference.rao.headings) < 1.0e-12, 1);

    figure_handle = figure('Color', 'w', 'Position', [60, 60, 1500, 900], ...
        'Name', '双球多体水动力与WAMIT对比');
    layout = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    plot_pair(layout, omega, squeeze(results.added_mass(1, 1, :)), ...
        squeeze(reference.radiation.added_mass(1, 1, reference_index)), ...
        squeeze(results.added_mass(1, 7, :)), ...
        squeeze(reference.radiation.added_mass(1, 7, reference_index)), ...
        '纵荡附加质量自项/互项', 'A_{11}, A_{1,7}（kg）');
    plot_pair(layout, omega, squeeze(results.added_mass(3, 3, :)), ...
        squeeze(reference.radiation.added_mass(3, 3, reference_index)), ...
        squeeze(results.added_mass(3, 9, :)), ...
        squeeze(reference.radiation.added_mass(3, 9, reference_index)), ...
        '升沉附加质量自项/互项', 'A_{33}, A_{3,9}（kg）');
    plot_pair(layout, omega, squeeze(results.damping(1, 1, :)), ...
        squeeze(reference.radiation.damping(1, 1, reference_index)), ...
        squeeze(results.damping(1, 7, :)), ...
        squeeze(reference.radiation.damping(1, 7, reference_index)), ...
        '纵荡阻尼自项/互项', 'B_{11}, B_{1,7}（kg/s）');
    plot_pair(layout, omega, squeeze(results.damping(3, 3, :)), ...
        squeeze(reference.radiation.damping(3, 3, reference_index)), ...
        squeeze(results.damping(3, 9, :)), ...
        squeeze(reference.radiation.damping(3, 9, reference_index)), ...
        '升沉阻尼自项/互项', 'B_{33}, B_{3,9}（kg/s）');

    nexttile(layout);
    hold on;
    plot(omega, squeeze(results.rao.amplitude(3, heading_zero, :)), 'b-o', ...
        'DisplayName', 'CRESTU-球1升沉');
    plot(omega, squeeze(results.rao.amplitude(9, heading_zero, :)), 'r-s', ...
        'DisplayName', 'CRESTU-球2升沉');
    plot(omega, squeeze(reference.rao.amplitude(3, reference_heading_zero, reference_index)), 'b--', ...
        'DisplayName', 'WAMIT-球1升沉');
    plot(omega, squeeze(reference.rao.amplitude(9, reference_heading_zero, reference_index)), 'r--', ...
        'DisplayName', 'WAMIT-球2升沉');
    format_axis('零度浪向两球升沉RAO', '|\Xi_3|, |\Xi_9|（m/m）');

    nexttile(layout);
    crestu_rao = squeeze(results.rao.amplitude(:, heading_zero, :)).';
    wamit_rao = squeeze(reference.rao.amplitude(:, reference_heading_zero, reference_index)).';
    dof_scale = max(max(abs(wamit_rao), [], 1), 1.0e-6);
    relative_error = 100 * abs(crestu_rao - wamit_rao) ./ dof_scale;
    imagesc(1:12, omega, relative_error);
    set(gca, 'YDir', 'normal');
    colorbar;
    grid on;
    xlabel('全局自由度编号', 'FontWeight', 'bold');
    ylabel('\omega（rad/s）', 'FontWeight', 'bold');
    title('12自由度RAO相对误差（%）', 'FontWeight', 'bold');

    title(layout, 'CRESTU双球相互作用与WAMIT TwoSphereIRR3对比', 'FontWeight', 'bold');
    output_file = fullfile(fileparts(summary_file), 'Plot_MultiBody_Results.png');
    export_chinese_figure(figure_handle, output_file);
    fprintf('>>> 已导出双球对比图: %s\n', output_file);
    if nargout == 0
        clear figure_handle
    end
end

function plot_pair(layout, omega, self_value, self_reference, cross_value, cross_reference, plot_title, y_label)
% PLOT_PAIR Plot one self coefficient and one cross-coupling coefficient.
%
% Syntax:
%   plot_pair(layout, omega, self_value, self_reference, cross_value, cross_reference, plot_title, y_label)
%
% Inputs:
%   layout          : [TiledChartLayout] Parent layout.
%   omega           : [1 x Nf] Angular frequencies, in rad/s.
%   self_value      : [Nf x 1] CRESTU self coefficient.
%   self_reference  : [Nf x 1] WAMIT self coefficient.
%   cross_value     : [Nf x 1] CRESTU cross coefficient.
%   cross_reference : [Nf x 1] WAMIT cross coefficient.
%   plot_title      : [char|string] Chinese subplot title.
%   y_label         : [char|string] Chinese vertical-axis label.
%
% Outputs:
%   None; one tile is populated.
%
% Mathematical Reference:
%   Direct dimensional coefficient comparison.
    nexttile(layout);
    hold on;
    plot(omega, self_value, 'b-o', 'DisplayName', 'CRESTU自项');
    plot(omega, self_reference, 'b--', 'DisplayName', 'WAMIT自项');
    plot(omega, cross_value, 'r-s', 'DisplayName', 'CRESTU互项');
    plot(omega, cross_reference, 'r--', 'DisplayName', 'WAMIT互项');
    format_axis(plot_title, y_label);
end

function format_axis(plot_title, y_label)
% FORMAT_AXIS Apply consistent Chinese multi-body plot formatting.
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
    grid on;
    box on;
    xlabel('\omega（rad/s）', 'FontWeight', 'bold');
    ylabel(y_label, 'FontWeight', 'bold');
    title(plot_title, 'FontWeight', 'bold');
    legend('Location', 'best');
end
