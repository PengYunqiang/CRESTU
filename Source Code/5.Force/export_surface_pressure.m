function pressure_data = export_surface_pressure(filename,mesh,phi_total,omega,rho)
% EXPORT_SURFACE_PRESSURE Execute the documented export_surface_pressure operation.
%
% Syntax:
%   pressure_data = export_surface_pressure(filename,mesh,phi_total,omega,rho)
%
% Inputs:
%   filename        : [char|string] Input or output file path.
%   mesh            : [struct] Boundary mesh with geometry expressed in SI units.
%   phi_total       : [documented value] Input required by the implemented function contract.
%   omega           : [scalar] Angular frequency, in rad/s.
%   rho             : [scalar] Water density, in kg/m^3.
%
% Outputs:
%   pressure_data   : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
%EXPORT_SURFACE_PRESSURE Compute and save complex panel pressure (.mat/.vtk).
% Pressure follows exp(i*omega*t): p=-i*omega*rho*phi_total.
    if ~isstruct(mesh)||~isfield(mesh,'centers')||~isfield(mesh,'vertices')||~isfield(mesh,'n_panels')
        error('CRESTU:PressureMesh','mesh requires centers, vertices, and n_panels.');
    end
    if size(phi_total,1)~=mesh.n_panels
        error('CRESTU:PressureShape','Potential rows must equal body panel count.');
    end
    pressure=-1i*omega*rho*phi_total;
    pressure_data=struct('omega',omega,'potential',phi_total,'pressure',pressure, ...
        'pressure_real',real(pressure),'pressure_imag',imag(pressure), ...
        'pressure_amplitude',abs(pressure),'pressure_phase_deg',angle(pressure)*180/pi, ...
        'centers',mesh.centers,'vertices',mesh.vertices);
    if nargin<1||isempty(filename), return; end
    [folder,~,ext]=fileparts(filename);
    if ~isempty(folder)&&~exist(folder,'dir'), mkdir(folder); end
    switch lower(ext)
        case '.mat'
            save(filename,'pressure_data','-v7');
        case '.vtk'
            if size(phi_total,2)~=1
                error('CRESTU:VTKWaveCase','VTK export accepts one potential/heading column at a time.');
            end
            write_vtk(filename,mesh,pressure_data);
        otherwise
            error('CRESTU:PressureFormat','Pressure output extension must be .mat or .vtk.');
    end
    fprintf('>>> Surface pressure exported: %s\n',filename);
end

function write_vtk(filename,mesh,data)
% WRITE_VTK Execute the documented write_vtk operation.
%
% Syntax:
%   write_vtk(filename,mesh,data)
%
% Inputs:
%   filename        : [char|string] Input or output file path.
%   mesh            : [struct] Boundary mesh with geometry expressed in SI units.
%   data            : [documented value] Input required by the implemented function contract.
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
    fid=fopen(filename,'wt');
    if fid<0, error('CRESTU:PressureOpen','Cannot create %s.',filename); end
    cleanup=onCleanup(@()fclose(fid)); np=mesh.n_panels; vertices=reshape(mesh.vertices,np,4,3);
    fprintf(fid,'# vtk DataFile Version 3.0\nCRESTU complex dynamic pressure\nASCII\n');
    fprintf(fid,'DATASET UNSTRUCTURED_GRID\nPOINTS %d double\n',4*np);
    for p=1:np
        for v=1:4, fprintf(fid,'%.15g %.15g %.15g\n',vertices(p,v,1),vertices(p,v,2),vertices(p,v,3)); end
    end
    fprintf(fid,'CELLS %d %d\n',np,5*np);
    for p=1:np, q=4*(p-1); fprintf(fid,'4 %d %d %d %d\n',q,q+1,q+2,q+3); end
    fprintf(fid,'CELL_TYPES %d\n',np); fprintf(fid,'%d\n',repmat(9,np,1));
    fprintf(fid,'CELL_DATA %d\n',np);
    write_scalar(fid,'pressure_real',data.pressure_real);
    write_scalar(fid,'pressure_imag',data.pressure_imag);
    write_scalar(fid,'pressure_amplitude',data.pressure_amplitude);
    write_scalar(fid,'pressure_phase_deg',data.pressure_phase_deg);
end

function write_scalar(fid,name,values)
% WRITE_SCALAR Execute the documented write_scalar operation.
%
% Syntax:
%   write_scalar(fid,name,values)
%
% Inputs:
%   fid             : [documented value] Input required by the implemented function contract.
%   name            : [integer scalar or array] Discrete count or index required by the algorithm.
%   values          : [documented value] Input required by the implemented function contract.
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
    fprintf(fid,'SCALARS %s double 1\nLOOKUP_TABLE default\n',name); fprintf(fid,'%.15g\n',values);
end
