function reference = read_wamit_first_order(file_one, rho)
% READ_WAMIT_FIRST_ORDER Read wamit first order for the CRESTU hydrodynamic workflow.
%
% Syntax:
%   reference = read_wamit_first_order(file_one, rho)
%
% Description:
%   Computes hydrodynamic coefficients, loads, or rigid-body response.
%   Results retain the CRESTU global 6-DOF order.
%
% Inputs:
%   file_one           - [character vector or string scalar] WAMIT first-order output-file path.
%   rho                - [scalar] Fluid density, [kg/m^3].
%
% Outputs:
%   reference          - [struct] Imported reference hydrodynamic data in documented SI normalization.
%
% Governing Equations / Theory:
%   Linear unsteady Bernoulli pressure, generalized surface integration, hydrostatics, radiation energy, or the frequency-domain rigid-body equation as applicable.
%
% References:
%   - Newman, J. N. (1977), Marine Hydrodynamics; Faltinsen, O. M. (1990), Sea Loads on Ships and Offshore Structures.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    raw = readmatrix(file_one,'FileType','text');
    if size(raw, 2) < 5
        error('CRESTU:WamitFormat','Expected five columns in %s.', file_one);
    end
    periods = unique(raw(:, 1),'stable');
    nf = numel(periods);
    ndof = max(max(raw(:, 2:3)));
    A = zeros(ndof, ndof, nf);
    B = zeros(ndof, ndof, nf);
    omegas = 2 * pi ./ periods(:).';
    for rowIndex = 1:size(raw, 1)
        k = find(periods == raw(rowIndex, 1), 1);
        i = raw(rowIndex, 2);
        j = raw(rowIndex, 3);
        A(i, j, k) = rho * raw(rowIndex, 4);
        B(i, j, k) = rho * omegas(k) * raw(rowIndex, 5);
    end
    reference = struct('file', file_one,'periods', periods(:).','omegas', omegas, ...
'added_mass', A,'damping', B,'rho', rho);
end
