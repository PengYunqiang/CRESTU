function fig = plot_surface_pressure(mesh,pressure,component)
% PLOT_SURFACE_PRESSURE Execute the documented plot_surface_pressure operation.
%
% Syntax:
%   fig = plot_surface_pressure(mesh,pressure,component)
%
% Inputs:
%   mesh            : [struct] Boundary mesh with geometry expressed in SI units.
%   pressure        : [documented value] Input required by the implemented function contract.
%   component       : [documented value] Input required by the implemented function contract.
%
% Outputs:
%   fig             : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
%PLOT_SURFACE_PRESSURE Plot real, imaginary, amplitude, or phase pressure.
    if nargin<3||isempty(component), component='amplitude'; end
    switch lower(component)
        case 'real', values=real(pressure);
        case 'imag', values=imag(pressure);
        case {'amplitude','abs'}, values=abs(pressure);
        case {'phase','phase_deg'}, values=angle(pressure)*180/pi;
        otherwise, error('CRESTU:PressureComponent','Unknown pressure component: %s',component);
    end
    if size(values,2)~=1||numel(values)~=mesh.n_panels
        error('CRESTU:PressurePlotShape','Plot one pressure column with one value per panel.');
    end
    vertices=reshape(mesh.vertices,mesh.n_panels,4,3);
    X=reshape(vertices(:,:,1),mesh.n_panels,4).';
    Y=reshape(vertices(:,:,2),mesh.n_panels,4).';
    Z=reshape(vertices(:,:,3),mesh.n_panels,4).';
    fig=figure('Color','w','Name','CRESTU surface pressure');
    patch(X,Y,Z,values(:).','FaceColor','flat','EdgeColor',[.25,.25,.25]);
    axis equal; grid on; view(3); colorbar; xlabel('x'); ylabel('y'); zlabel('z');
    title(sprintf('Surface pressure: %s',component),'Interpreter','none');
end
