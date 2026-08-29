function massMatrix = assemble_mass_matrix(mass_props)
% ASSEMBLE_MASS_MATRIX Assemble the global rigid-body mass and inertia matrix.
%
% Syntax:
%   massMatrix = assemble_mass_matrix(mass_props)
%
% Description:
%   Computes hydrodynamic coefficients, loads, or rigid-body response.
%   Results retain the CRESTU global 6-DOF order.
%
% Inputs:
%   mass_props         - [struct array] Body masses, centers of gravity, and inertia tensors in SI units.
%
% Outputs:
%   massMatrix         - [6Nb x 6Nb] Global rigid-body mass and inertia matrix, [kg], [kg m], and [kg m^2].
%
% Governing Equations / Theory:
%   Linear unsteady Bernoulli pressure, generalized surface integration, hydrostatics, radiation energy, or the frequency-domain rigid-body equation as applicable.
%
% References:
%   - Newman, J. N. (1977), Marine Hydrodynamics; Faltinsen, O. M. (1990), Sea Loads on Ships and Offshore Structures.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    bodyCount = numel(mass_props);
    assert(bodyCount > 0, 'CRESTU:MassProperties', ...
        'At least one body mass-property record is required.');

    %% Stage 2: Assemble body mass and inertia blocks

    % A five-module model has 30 global DOFs. Module m uses
    % 6*(m-1)+(1:6). Local DOF 1-6: Surge, Sway, Heave, Roll, Pitch, Yaw.
    massMatrix = zeros(6 * bodyCount);

    for bodyIndex = 1:bodyCount
        globalDofRows = 6 * (bodyIndex - 1) + (1:6);
        bodyMass = mass_props(bodyIndex).mass; % [kg]
        bodyInertia = mass_props(bodyIndex).inertia; % [kg*m^2]

        if ~isequal(size(bodyInertia), [3, 3]) || bodyMass <= 0 || ...
                any(eig((bodyInertia + bodyInertia.') / 2) <= 0)
            error('CRESTU:MassProperties', ...
                'Invalid mass properties for body %d.', bodyIndex);
        end

        massMatrix(globalDofRows, globalDofRows) = ...
            blkdiag(bodyMass * eye(3), (bodyInertia + bodyInertia.') / 2);
    end
end
