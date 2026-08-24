function export_chinese_figure(figure_handle, filename)
% EXPORT_CHINESE_FIGURE Export a publication-quality Chinese-labelled figure.
%
% Syntax:
%   export_chinese_figure(figure_handle, filename)
%
% Inputs:
%   figure_handle : [Figure] MATLAB figure to export.
%   filename      : [char|string] PNG output path.
%
% Outputs:
%   None; a 300-DPI PNG file is written.
%
% Mathematical Reference:
%   Visualization utility; no mathematical model is used.
    output_folder = fileparts(filename);
    if ~isempty(output_folder) && ~isfolder(output_folder)
        mkdir(output_folder);
    end
    set(findall(figure_handle, '-property', 'FontName'), 'FontName', 'Microsoft YaHei');
    set(findall(figure_handle, '-property', 'FontSize'), 'FontSize', 10);
    axes_handles = findall(figure_handle, 'Type', 'axes');
    for axes_index = 1:numel(axes_handles)
        disableDefaultInteractivity(axes_handles(axes_index));
        axes_handles(axes_index).Toolbar.Visible = 'off';
    end
    exportgraphics(figure_handle, filename, 'Resolution', 300);
end
