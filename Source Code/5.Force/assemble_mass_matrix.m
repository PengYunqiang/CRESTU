function mass_matrix = assemble_mass_matrix(mass_props)
% ASSEMBLE_MASS_MATRIX Assemble the global rigid-body mass and inertia matrix.
%
% Syntax:
%   mass_matrix = assemble_mass_matrix(mass_props)
%
% Description:
%   The routine converts first-order potential-flow or rigid-body data into generalized hydrodynamic coefficients, loads, restoring terms, or motions. Translational and rotational quantities retain the global 6-DOF ordering used by CRESTU.
%
% Inputs:
%   mass_props         - [struct array] Body masses, centers of gravity, and inertia tensors in SI units.
%
% Outputs:
%   mass_matrix        - [6Nb x 6Nb] Global rigid-body mass and inertia matrix, [kg], [kg m], and [kg m^2].
%
% Governing Equations / Theory:
%   Linear unsteady Bernoulli pressure, generalized surface integration, hydrostatics, radiation energy, or the frequency-domain rigid-body equation as applicable.
%
% References:
%   - Newman, J. N. (1977), Marine Hydrodynamics; Faltinsen, O. M. (1990), Sea Loads on Ships and Offshore Structures.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% --- 1. Validate Inputs and Initialize the Algorithm ---

    n = numel(mass_props); mass_matrix = zeros(6 * n);
    for b = 1:n
        idx = (b - 1) * 6 + (1:6); m = mass_props(b).mass; inertia = mass_props(b).inertia;
        if ~isequal(size(inertia), [3, 3]) || m <= 0 || any(eig((inertia + inertia.') / 2) <= 0)
            error('CRESTU:MassProperties', 'Invalid mass properties for body %d.', b);
        end
        mass_matrix(idx, idx) = blkdiag(m * eye(3), (inertia + inertia.') / 2);
    end
end
