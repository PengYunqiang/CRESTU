function parity = get_mode_parities(n_bodies, isx, isy)
% GET_MODE_PARITIES Return the reflection parity of every rigid-body radiation mode.
%
% Syntax:
%   parity = get_mode_parities(n_bodies, isx, isy)
%
% Description:
%   Implements linear Rankine potential-flow operations.
%   Symmetry and phase follow exp(i*omega*t).
%
% Inputs:
%   n_bodies           - [scalar] Number of hydrodynamically coupled bodies, dimensionless.
%   isx                - [logical scalar] Symmetry-reduction flag for the x = 0 plane.
%   isy                - [logical scalar] Symmetry-reduction flag for the y = 0 plane.
%
% Outputs:
%   parity             - [1 x 2] Reflection parity signs for the x = 0 and y = 0 planes, dimensionless.
%
% Governing Equations / Theory:
%   Green third identity, the Rankine kernel 1/r, linearized free-surface theory, and reflection symmetry as applicable.
%
% References:
%   - Newman, J. N. (1977), Marine Hydrodynamics; Hess and Smith (1964); project boundary-condition specification.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    px = [-1, 1, 1, 1, -1, -1];
    py = [1, -1, 1, -1, 1, -1];
    parity = [repmat(px(:), n_bodies, 1), repmat(py(:), n_bodies, 1)];
    if ~isx
        parity(:, 1) = 1;
    end
    if ~isy
        parity(:, 2) = 1;
    end
end
