function comparison = compare_wamit_results(results, wamit_reference)
% COMPARE_WAMIT_RESULTS Compare wamit results for the CRESTU hydrodynamic workflow.
%
% Syntax:
%   comparison = compare_wamit_results(results, wamit_reference)
%
% Description:
%   The routine converts first-order potential-flow or rigid-body data into generalized hydrodynamic coefficients, loads, restoring terms, or motions. Translational and rotational quantities retain the global 6-DOF ordering used by CRESTU.
%
% Inputs:
%   results            - [struct] CRESTU hydrodynamic results in the project SI normalization.
%   wamit_reference    - [struct] WAMIT reference results converted to the project SI normalization.
%
% Outputs:
%   comparison         - [struct] Value, reference, and relative-error comparison fields.
%
% Governing Equations / Theory:
%   Linear unsteady Bernoulli pressure, generalized surface integration, hydrostatics, radiation energy, or the frequency-domain rigid-body equation as applicable.
%
% References:
%   - Newman, J. N. (1977), Marine Hydrodynamics; Faltinsen, O. M. (1990), Sea Loads on Ships and Offshore Structures.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% --- 1. Validate Inputs and Initialize the Algorithm ---

    nf = numel(results.omegas);
    template = struct('omega', 0, 'reference_omega', 0, 'A33', 0, 'A33_reference', 0, ...
        'A33_relative_error', 0, 'B33', 0, 'B33_reference', 0, 'B33_relative_error', 0, ...
        'A_symmetry_error', 0, 'B_symmetry_error', 0);
    comparison = repmat(template, nf, 1);
    for k = 1:nf
        [~, q] = min(abs(wamit_reference.omegas - results.omegas(k)));
        A = results.added_mass(:, :, k); B = results.damping(:, :, k);
        Ar = wamit_reference.added_mass(:, :, q); Br = wamit_reference.damping(:, :, q);
        comparison(k) = struct('omega', results.omegas(k), 'reference_omega', wamit_reference.omegas(q), ...
            'A33', A(3, 3), 'A33_reference', Ar(3, 3), 'A33_relative_error', relative_error(A(3, 3), Ar(3, 3)), ...
            'B33', B(3, 3), 'B33_reference', Br(3, 3), 'B33_relative_error', relative_error(B(3, 3), Br(3, 3)), ...
            'A_symmetry_error', norm(A - A.', 'fro') / max(norm(A, 'fro'), eps), ...
            'B_symmetry_error', norm(B - B.', 'fro') / max(norm(B, 'fro'), eps));
    end
end

function e = relative_error(value, reference)
% RELATIVE_ERROR Evaluate a scale-aware relative error.
%
% Syntax:
%   e = relative_error(value, reference)
%
% Description:
%   The routine converts first-order potential-flow or rigid-body data into generalized hydrodynamic coefficients, loads, restoring terms, or motions. Translational and rotational quantities retain the global 6-DOF ordering used by CRESTU.
%
% Inputs:
%   value              - [numeric scalar or array] Computed quantity to compare with a reference value, in matching SI units.
%   reference          - [struct or numeric array] Reference solution and normalization metadata in stated SI units.
%
% Outputs:
%   e                  - [numeric array] Scale-aware relative error, dimensionless.
%
% Governing Equations / Theory:
%   Linear unsteady Bernoulli pressure, generalized surface integration, hydrostatics, radiation energy, or the frequency-domain rigid-body equation as applicable.
%
% References:
%   - Newman, J. N. (1977), Marine Hydrodynamics; Faltinsen, O. M. (1990), Sea Loads on Ships and Offshore Structures.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% --- 1. Validate Inputs and Initialize the Algorithm ---

    e = abs(value - reference) / max(abs(reference), eps);
end
