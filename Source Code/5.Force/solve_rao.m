function rao = solve_rao(omegas, added_mass, damping, hydrostatic, excitation, cfg)
% SOLVE_RAO Solve the coupled frequency-domain rigid-body equations for response amplitude operators.
%
% Syntax:
%   rao = solve_rao(omegas, added_mass, damping, hydrostatic, excitation, cfg)
%
% Description:
%   The routine converts first-order potential-flow or rigid-body data into generalized hydrodynamic coefficients, loads, restoring terms, or motions. Translational and rotational quantities retain the global 6-DOF ordering used by CRESTU.
%
% Inputs:
%   omegas             - [1 x Nf] Angular-frequency samples, [rad/s].
%   added_mass         - [Ndof x Ndof x Nf] Added-mass coefficients in consistent translational and rotational SI units.
%   damping            - [Ndof x Ndof x Nf] Radiation-damping coefficients in consistent translational and rotational SI units.
%   hydrostatic        - [Ndof x Ndof] Linear hydrostatic restoring matrix in consistent SI units.
%   excitation         - [Ndof x Nh x Nf] Complex first-order wave-excitation loads, [N] and [N m].
%   cfg                - [struct] Validated CRESTU configuration containing SI-valued physical and numerical parameters.
%
% Outputs:
%   rao                - [Ndof x Nh x Nf] Complex motion response per unit wave amplitude, [m/m] and [rad/m].
%
% Governing Equations / Theory:
%   Linear unsteady Bernoulli pressure, generalized surface integration, hydrostatics, radiation energy, or the frequency-domain rigid-body equation as applicable.
%
% References:
%   - Newman, J. N. (1977), Marine Hydrodynamics; Faltinsen, O. M. (1990), Sea Loads on Ships and Offshore Structures.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% --- 1. Validate Inputs and Initialize the Algorithm ---

    omegas = reshape(omegas, 1, []); nf = numel(omegas); ndof = 6 * cfg.n_bodies;
    added_mass = normalize_coeff(added_mass, ndof, nf, 'added mass');
    damping = normalize_coeff(damping, ndof, nf, 'damping');
    if ~isequal(size(hydrostatic), [ndof, ndof])
        error('CRESTU:HydrostaticShape', 'Hydrostatic matrix must be %d-by-%d.', ndof, ndof);
    end
    if ismatrix(excitation) && nf == 1, excitation = reshape(excitation, ndof, size(excitation, 2), 1); end
    if size(excitation, 1) ~= ndof || size(excitation, 3) ~= nf
        error('CRESTU:ExcitationShape', 'Excitation must be ndof-by-nheading-by-nfrequency.');
    end
    nh = size(excitation, 2); M = assemble_mass_matrix(cfg.mass_props);
    xi = complex(zeros(ndof, nh, nf)); dynamic_matrix = complex(zeros(ndof, ndof, nf)); rconds = zeros(nf, 1);
    for k = 1:nf
        w = omegas(k); Z = -w^2 * (M + added_mass(:, :, k)) + 1i * w * damping(:, :, k) + hydrostatic;
        dynamic_matrix(:, :, k) = Z; rconds(k) = rcond(Z);
        if rconds(k) < 1e-12, warning('CRESTU:IllConditionedRAO', 'RAO matrix at omega=%g has rcond=%g.', w, rconds(k)); end
        xi(:, :, k) = Z \ excitation(:, :, k);
    end
    rao = struct('complex', xi, 'amplitude', abs(xi), 'phase_deg', angle(xi) * 180 / pi, ...
        'mass_matrix', M, 'dynamic_matrix', dynamic_matrix, 'rcond', rconds, ...
        'omegas', omegas, 'headings', cfg.wave.headings);
end

function out = normalize_coeff(in, n, nf, label)
% NORMALIZE_COEFF Normalize frequency-dependent coefficient arrays to the required shape.
%
% Syntax:
%   out = normalize_coeff(in, n, nf, label)
%
% Description:
%   The routine converts first-order potential-flow or rigid-body data into generalized hydrodynamic coefficients, loads, restoring terms, or motions. Translational and rotational quantities retain the global 6-DOF ordering used by CRESTU.
%
% Inputs:
%   in                 - [character vector, string scalar, or numeric value] Value to normalize.
%   n                  - [scalar] Number of coupled degrees of freedom or array rows, dimensionless.
%   nf                 - [scalar] Number of frequency samples, dimensionless.
%   label              - [character vector or string scalar] Diagnostic field label.
%
% Outputs:
%   out                - [value] Normalized or selected result with type and physical units preserved from the input.
%
% Governing Equations / Theory:
%   Linear unsteady Bernoulli pressure, generalized surface integration, hydrostatics, radiation energy, or the frequency-domain rigid-body equation as applicable.
%
% References:
%   - Newman, J. N. (1977), Marine Hydrodynamics; Faltinsen, O. M. (1990), Sea Loads on Ships and Offshore Structures.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% --- 1. Validate Inputs and Initialize the Algorithm ---

    if iscell(in), out = zeros(n, n, nf); for k = 1:nf, out(:, :, k) = in{k}; end
    else, out = in; end
    if nf == 1 && ismatrix(out), out = reshape(out, n, n, 1); end
    if size(out, 1) ~= n || size(out, 2) ~= n || size(out, 3) ~= nf
        error('CRESTU:CoefficientShape', '%s array has the wrong shape.', label);
    end
end
