function write_bmf(filename, mesh)
% WRITE_BMF Write bmf for the CRESTU hydrodynamic workflow.
%
% Syntax:
%   write_bmf(filename, mesh)
%
% Description:
%   Prepares boundary-panel geometry for the Rankine solver.
%   Global coordinates and panel-normal signs are preserved.
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

%% 阶段 1: 以短时重试方式打开输出文件

    [fid, openMessage] = open_bmf_with_retry(filename, 10, 0.05); % [s]

    if fid == -1
        error('CRESTU:BmfWriteOpen', ...
            'Unable to create the BMF file: %s (%s)', filename, openMessage);
    end

    fileCleanup = onCleanup(@() fclose(fid));

%% 阶段 2: 写入 BMF 头部和逐面元顶点

    unique_types = unique(mesh.panel_type);
    if length(unique_types) == 1
        main_type = unique_types(1);
    else
        main_type = 0;
    end

    fprintf(fid,'%s\n', mesh.header);
    fprintf(fid,'%16.6E %5d\n', mesh.ulen, main_type);
    fprintf(fid,'%5d %5d\n', mesh.isx, mesh.isy);
    fprintf(fid,'%8d\n', mesh.n_panels);

    for k = 1:mesh.n_panels
        for vertexIndex = 1:4
            fprintf(fid,'%18.9E %18.9E %18.9E\n', ...
                mesh.vertices(k, vertexIndex, 1), mesh.vertices(k, vertexIndex, 2), mesh.vertices(k, vertexIndex, 3));
        end
    end
    clear fileCleanup

    fprintf('[OK] BMF mesh exported: %s (NPAN=%d, Type=%d)\n', filename, mesh.n_panels, main_type);
end

function [fileIdentifier, openMessage] = open_bmf_with_retry( ...
        filename, maximumAttemptCount, retryDelaySeconds)
% OPEN_BMF_WITH_RETRY Tolerate short-lived Windows file-indexing locks.

    fileIdentifier = -1;
    openMessage = '';

    for attemptIndex = 1:maximumAttemptCount
        [fileIdentifier, openMessage] = fopen(filename, 'w');

        if fileIdentifier >= 0
            return
        end

        if attemptIndex < maximumAttemptCount
            fprintf('[WARN] BMF write retry %d/%d | %s\n', ...
                attemptIndex, maximumAttemptCount, filename);
            pause(retryDelaySeconds); % [s]
        end
    end
end
