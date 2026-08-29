function expanded = expand_scalar_by_symmetry(values, parity, isx, isy)
% EXPAND_SCALAR_BY_SYMMETRY Expand scalar by symmetry for the CRESTU hydrodynamic workflow.
%
% Syntax:
%   expanded = expand_scalar_by_symmetry(values, parity, isx, isy)
%
% Description:
%   Computes or processes second-order mean wave-drift quantities.
%   Complex averages follow exp(i*omega*t) and the 6-DOF order.
%
% Inputs:
%   values             - [numeric array] Samples to be transformed or interpolated; units are preserved.
%   parity             - [1 x 2] Reflection parity signs for the x = 0 and y = 0 planes, dimensionless.
%   isx                - [logical scalar] Symmetry-reduction flag for the x = 0 plane.
%   isy                - [logical scalar] Symmetry-reduction flag for the y = 0 plane.
%
% Outputs:
%   expanded           - [numeric array] Full-domain scalar samples with the same units as the reduced input.
%
% Governing Equations / Theory:
%   Pinkster near-field direct pressure integration, Maruo-Newman far-field momentum balance, or supporting surface reconstruction as applicable.
%
% References:
%   - Pinkster, J. A. (1980), Low Frequency Second Order Wave Exciting Forces on Floating Structures; Newman, J. N. (1974), second-order slowly varying forces.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    if isempty(values)
        expanded = values;
        return;
    end
    parity = reshape(parity, size(values, 2), 2);
    flags = [0, 0];
    if isx
        flags = [flags;1, 0];
    end
    if isy
        flags = [flags;0, 1];
    end
    if isx && isy
        flags = [0, 0;1, 0;0, 1;1, 1];
    end
    expanded = complex(zeros(size(values, 1) * size(flags, 1), size(values, 2)));
    for imageIndex = 1:size(flags, 1)
        rows = (imageIndex - 1) * size(values, 1) + (1:size(values, 1));
        weights = (parity(:, 1) .^ flags(imageIndex, 1)) .* (parity(:, 2) .^ flags(imageIndex, 2));
        expanded(rows, :) = values .* weights.';
    end
end
