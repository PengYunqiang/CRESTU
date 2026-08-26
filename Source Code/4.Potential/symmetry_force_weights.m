function weights = symmetry_force_weights(output_parity, solution_parity, isx, isy)
% SYMMETRY_FORCE_WEIGHTS Construct force-recovery weights for symmetry-reduced solutions.
%
% Syntax:
%   weights = symmetry_force_weights(output_parity, solution_parity, isx, isy)
%
% Description:
%   The routine implements a component of the linear Rankine boundary-element formulation for incompressible, irrotational gravity-wave flow. Geometry, reflection parity, free-surface impedance, and complex phase follow the project convention exp(i*omega*t).
%
% Inputs:
%   output_parity      - [1 x 2] Reflection parity of the recovered load component, dimensionless.
%   solution_parity    - [1 x 2] Reflection parity of the solved potential component, dimensionless.
%   isx                - [logical scalar] Symmetry-reduction flag for the x = 0 plane.
%   isy                - [logical scalar] Symmetry-reduction flag for the y = 0 plane.
%
% Outputs:
%   weights            - [numeric array] Symmetry, interpolation, or quadrature weights, dimensionless.
%
% Governing Equations / Theory:
%   Green third identity, the Rankine kernel 1/r, linearized free-surface theory, and reflection symmetry as applicable.
%
% References:
%   - Newman, J. N. (1977), Marine Hydrodynamics; Hess and Smith (1964); project boundary-condition specification.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% --- 1. Validate Inputs and Initialize the Algorithm ---

    weights = ones(size(output_parity, 1), size(solution_parity, 1));
    if isx, weights = weights .* (1 + output_parity(:, 1) * solution_parity(:, 1).'); end
    if isy, weights = weights .* (1 + output_parity(:, 2) * solution_parity(:, 2).'); end
end
