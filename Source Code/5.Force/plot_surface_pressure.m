function fig = plot_surface_pressure(mesh, pressure, component)
% PLOT_SURFACE_PRESSURE Plot surface pressure for the CRESTU hydrodynamic workflow.
%
% Syntax:
%   fig = plot_surface_pressure(mesh, pressure, component)
%
% Description:
%   The routine converts first-order potential-flow or rigid-body data into generalized hydrodynamic coefficients, loads, restoring terms, or motions. Translational and rotational quantities retain the global 6-DOF ordering used by CRESTU.
%
% Inputs:
%   mesh               - [struct] Boundary-panel mesh with Cartesian geometry in SI units.
%   pressure           - [N x K] Complex panel pressure or a selected pressure component, [Pa].
%   component          - [character vector or string scalar] Pressure component selector.
%
% Outputs:
%   fig                - [graphics handle] MATLAB figure containing the requested visualization.
%
% Governing Equations / Theory:
%   Linear unsteady Bernoulli pressure, generalized surface integration, hydrostatics, radiation energy, or the frequency-domain rigid-body equation as applicable.
%
% References:
%   - Newman, J. N. (1977), Marine Hydrodynamics; Faltinsen, O. M. (1990), Sea Loads on Ships and Offshore Structures.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% --- 1. Validate Inputs and Initialize the Algorithm ---

    if nargin < 3 || isempty(component), component = 'amplitude'; end
    switch lower(component)
        case 'real', values = real(pressure);
        case 'imag', values = imag(pressure);
        case {'amplitude', 'abs'}, values = abs(pressure);
        case {'phase', 'phase_deg'}, values = angle(pressure) * 180 / pi;
        otherwise, error('CRESTU:PressureComponent', 'Unknown pressure component: %s', component);
    end
    if size(values, 2) ~= 1 || numel(values) ~= mesh.n_panels
        error('CRESTU:PressurePlotShape', 'Plot one pressure column with one value per panel.');
    end
    vertices = reshape(mesh.vertices, mesh.n_panels, 4, 3);
    X = reshape(vertices(:, :, 1), mesh.n_panels, 4).';
    Y = reshape(vertices(:, :, 2), mesh.n_panels, 4).';
    Z = reshape(vertices(:, :, 3), mesh.n_panels, 4).';
    fig = figure('Color', 'w', 'Name', 'CRESTU surface pressure');
    patch(X, Y, Z, values(:).', 'FaceColor', 'flat', 'EdgeColor', [.25, .25, .25]);
    axis equal; grid on; view(3); colorbar; xlabel('x'); ylabel('y'); zlabel('z');
    title(sprintf('Surface pressure: %s', component), 'Interpreter', 'none');
end
