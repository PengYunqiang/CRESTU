function reference = read_wamit_rao(filename)
% READ_WAMIT_RAO Read wamit rao for the CRESTU hydrodynamic workflow.
%
% Syntax:
%   reference = read_wamit_rao(filename)
%
% Description:
%   Computes hydrodynamic coefficients, loads, or rigid-body response.
%   Results retain the CRESTU global 6-DOF order.
%
% Inputs:
%   filename           - [character vector or string scalar] Input or output file path.
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

    raw = readmatrix(filename,'FileType','text');
    if size(raw, 2) < 7
        error('CRESTU:WamitFormat','Expected seven columns in %s.', filename);
    end
    periods = unique(raw(:, 1),'stable');
    headings = unique(raw(:, 2),'stable');
    nf = numel(periods);
    nh = numel(headings);
    ndof = max(raw(:, 3));
    value = complex(zeros(ndof, nh, nf));
    for rowIndex = 1:size(raw, 1)
        k = find(periods == raw(rowIndex, 1), 1);
        headingIndex = find(headings == raw(rowIndex, 2), 1);
        mode = raw(rowIndex, 3);
        value(mode, headingIndex, k) = complex(raw(rowIndex, 6), raw(rowIndex, 7));
    end
    reference = struct('file', filename,'periods', periods(:).','omegas', 2 * pi ./ periods(:).', ...
'headings', headings(:).','complex', value,'amplitude', abs(value),'phase_deg', angle(value) * 180 / pi);
end
