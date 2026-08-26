function write_bmf(filename, mesh)
% WRITE_BMF Write bmf for the CRESTU hydrodynamic workflow.
%
% Syntax:
%   write_bmf(filename, mesh)
%
% Description:
%   The routine constructs, transforms, validates, or visualizes boundary-panel geometry used by the Rankine solver. Coordinates are expressed in the global Cartesian frame and panel orientation is preserved so that normals remain consistent with boundary-integral signs.
%
% Inputs:
%   filename           - [character vector or string scalar] Input or output file path.
%   mesh               - [struct] Boundary-panel mesh with Cartesian geometry in SI units.
%
% Outputs:
%   None.
%
% Governing Equations / Theory:
%   Planar panel geometry, polygon moments, reflection transformations, and mesh-topology relations as applicable.
%
% References:
%   - Hess and Smith (1964); CRESTU BMF mesh-format and geometry conventions.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% --- 1. Validate Inputs and Initialize the Algorithm ---

    fid = fopen(filename, 'w');
    if fid == -1
        error('Unable to create the BMF file: %s', filename);
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

    fprintf('>>> BMF mesh exported successfully: %s (NPAN=%d, Type=%d)\n', filename, mesh.n_panels, main_type);
end
