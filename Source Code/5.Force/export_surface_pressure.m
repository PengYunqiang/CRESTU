function pressure_data = export_surface_pressure(filename, mesh, phi_total, omega, rho)
% EXPORT_SURFACE_PRESSURE Export surface pressure for the CRESTU hydrodynamic workflow.
%
% Syntax:
%   pressure_data = export_surface_pressure(filename, mesh, phi_total, omega, rho)
%
% Description:
%   The routine converts first-order potential-flow or rigid-body data into generalized hydrodynamic coefficients, loads, restoring terms, or motions. Translational and rotational quantities retain the global 6-DOF ordering used by CRESTU.
%
% Inputs:
%   filename           - [character vector or string scalar] Input or output file path.
%   mesh               - [struct] Boundary-panel mesh with Cartesian geometry in SI units.
%   phi_total          - [N x Nh] Complex total first-order potentials, [m^2/s].
%   omega              - [scalar] Angular frequency, [rad/s].
%   rho                - [scalar] Fluid density, [kg/m^3].
%
% Outputs:
%   pressure_data      - [struct] Exported complex pressure components and file metadata, [Pa].
%
% Governing Equations / Theory:
%   Linear unsteady Bernoulli pressure, generalized surface integration, hydrostatics, radiation energy, or the frequency-domain rigid-body equation as applicable.
%
% References:
%   - Newman, J. N. (1977), Marine Hydrodynamics; Faltinsen, O. M. (1990), Sea Loads on Ships and Offshore Structures.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% --- 1. Validate Inputs and Initialize the Algorithm ---

    if ~isstruct(mesh) || ~isfield(mesh, 'centers') || ~isfield(mesh, 'vertices') || ~isfield(mesh, 'n_panels')
        error('CRESTU:PressureMesh', 'mesh requires centers, vertices, and n_panels.');
    end
    if size(phi_total, 1) ~= mesh.n_panels
        error('CRESTU:PressureShape', 'Potential rows must equal body panel count.');
    end
    pressure = -1i * omega * rho * phi_total;
    pressure_data = struct('omega', omega, 'potential', phi_total, 'pressure', pressure, ...
        'pressure_real', real(pressure), 'pressure_imag', imag(pressure), ...
        'pressure_amplitude', abs(pressure), 'pressure_phase_deg', angle(pressure) * 180 / pi, ...
        'centers', mesh.centers, 'vertices', mesh.vertices);
    if nargin < 1 || isempty(filename), return; end
    [folder, ~, ext] = fileparts(filename);
    if ~isempty(folder) && ~exist(folder, 'dir'), mkdir(folder); end
    switch lower(ext)
        case '.mat'
            save(filename, 'pressure_data', '-v7');
        case '.vtk'
            if size(phi_total, 2) ~= 1
                error('CRESTU:VTKWaveCase', 'VTK export accepts one potential/heading column at a time.');
            end
            write_vtk(filename, mesh, pressure_data);
        otherwise
            error('CRESTU:PressureFormat', 'Pressure output extension must be .mat or .vtk.');
    end
    fprintf('>>> Surface pressure exported: %s\n', filename);
end

function write_vtk(filename, mesh, data)
% WRITE_VTK Write panel geometry and scalar fields in legacy VTK format.
%
% Syntax:
%   write_vtk(filename, mesh, data)
%
% Description:
%   The routine converts first-order potential-flow or rigid-body data into generalized hydrodynamic coefficients, loads, restoring terms, or motions. Translational and rotational quantities retain the global 6-DOF ordering used by CRESTU.
%
% Inputs:
%   filename           - [character vector or string scalar] Input or output file path.
%   mesh               - [struct] Boundary-panel mesh with Cartesian geometry in SI units.
%   data               - [struct] Scalar panel fields to write, with units documented by each field.
%
% Outputs:
%   None.
%
% Governing Equations / Theory:
%   Linear unsteady Bernoulli pressure, generalized surface integration, hydrostatics, radiation energy, or the frequency-domain rigid-body equation as applicable.
%
% References:
%   - Newman, J. N. (1977), Marine Hydrodynamics; Faltinsen, O. M. (1990), Sea Loads on Ships and Offshore Structures.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% --- 1. Validate Inputs and Initialize the Algorithm ---

    fid = fopen(filename, 'wt');
    if fid < 0, error('CRESTU:PressureOpen', 'Cannot create %s.', filename); end
    cleanup = onCleanup(@()fclose(fid)); np = mesh.n_panels; vertices = reshape(mesh.vertices, np, 4, 3);
    fprintf(fid, '# vtk DataFile Version 3.0\nCRESTU complex dynamic pressure\nASCII\n');
    fprintf(fid, 'DATASET UNSTRUCTURED_GRID\nPOINTS %d double\n', 4 * np);
    for p = 1:np
        for v = 1:4, fprintf(fid, '%.15g %.15g %.15g\n', vertices(p, v, 1), vertices(p, v, 2), vertices(p, v, 3)); end
    end
    fprintf(fid, 'CELLS %d %d\n', np, 5 * np);
    for p = 1:np, q = 4 * (p - 1); fprintf(fid, '4 %d %d %d %d\n', q, q + 1, q + 2, q + 3); end
    fprintf(fid, 'CELL_TYPES %d\n', np); fprintf(fid, '%d\n', repmat(9, np, 1));
    fprintf(fid, 'CELL_DATA %d\n', np);
    write_scalar(fid, 'pressure_real', data.pressure_real);
    write_scalar(fid, 'pressure_imag', data.pressure_imag);
    write_scalar(fid, 'pressure_amplitude', data.pressure_amplitude);
    write_scalar(fid, 'pressure_phase_deg', data.pressure_phase_deg);
end

function write_scalar(fid, name, values)
% WRITE_SCALAR Write one scalar field to an open VTK file.
%
% Syntax:
%   write_scalar(fid, name, values)
%
% Description:
%   The routine converts first-order potential-flow or rigid-body data into generalized hydrodynamic coefficients, loads, restoring terms, or motions. Translational and rotational quantities retain the global 6-DOF ordering used by CRESTU.
%
% Inputs:
%   fid                - [scalar] Valid open MATLAB file identifier, dimensionless.
%   name               - [character vector or string scalar] Output scalar-field name.
%   values             - [numeric array] Samples to be transformed or interpolated; units are preserved.
%
% Outputs:
%   None.
%
% Governing Equations / Theory:
%   Linear unsteady Bernoulli pressure, generalized surface integration, hydrostatics, radiation energy, or the frequency-domain rigid-body equation as applicable.
%
% References:
%   - Newman, J. N. (1977), Marine Hydrodynamics; Faltinsen, O. M. (1990), Sea Loads on Ships and Offshore Structures.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% --- 1. Validate Inputs and Initialize the Algorithm ---

    fprintf(fid, 'SCALARS %s double 1\nLOOKUP_TABLE default\n', name); fprintf(fid, '%.15g\n', values);
end
