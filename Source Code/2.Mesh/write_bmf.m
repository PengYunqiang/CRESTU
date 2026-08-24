function write_bmf(filename, mesh)
% WRITE_BMF Execute the documented write_bmf operation.
%
% Syntax:
%   write_bmf(filename, mesh)
%
% Inputs:
%   filename        : [char|string] Input or output file path.
%   mesh            : [struct] Boundary mesh with geometry expressed in SI units.
%
% Outputs:
%   None; the function performs the documented file, plot, or validation action.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
%
%   Line 1: Header
%   Line 2: Length Scale,Panel Type
%   Line 3: ISX, ISY
%   Line 4: NPAN

    fid = fopen(filename, 'w');
    if fid == -1
        error('无法创建 BMF 文件: %s', filename);
    end

    unique_types = unique(mesh.panel_type);
    if length(unique_types) == 1
        main_type = unique_types(1);
    else
        main_type = 0;
    end

    fprintf(fid, '%s\n', mesh.header);
    fprintf(fid, '%16.6E %5d\n', mesh.ulen, main_type);
    fprintf(fid, '%5d %5d\n', mesh.isx, mesh.isy);
    fprintf(fid, '%8d\n', mesh.n_panels);

    for k = 1:mesh.n_panels
        for v = 1:4
            fprintf(fid, '%18.9E %18.9E %18.9E\n', ...
                mesh.vertices(k, v, 1), mesh.vertices(k, v, 2), mesh.vertices(k, v, 3));
        end
    end
    fclose(fid);

    fprintf('>>> BMF 网格已成功导出: %s (NPAN=%d, Type=%d)\n', filename, mesh.n_panels, main_type);
end
