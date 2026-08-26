function [X, factors] = solve_complex_system(A, B, factors)
% SOLVE_COMPLEX_SYSTEM Solve a reusable complex linear system with cached factorization.
%
% Syntax:
%   [X, factors] = solve_complex_system(A, B, factors)
%
% Description:
%   The routine implements a component of the linear Rankine boundary-element formulation for incompressible, irrotational gravity-wave flow. Geometry, reflection parity, free-surface impedance, and complex phase follow the project convention exp(i*omega*t).
%
% Inputs:
%   A                  - [numeric array] Matrix or geometric area defined by the governing formulation, in corresponding SI units.
%   B                  - [numeric array] Secondary matrix or geometric breadth defined by the function contract, in corresponding SI units.
%   factors            - [struct or empty] Reusable complex-system factorization data.
%
% Outputs:
%   X                  - [numeric array] Complex solution matrix with units inherited from the right-hand side.
%   factors            - [struct or empty] Reusable complex-system factorization data.
%
% Governing Equations / Theory:
%   Green third identity, the Rankine kernel 1/r, linearized free-surface theory, and reflection symmetry as applicable.
%
% References:
%   - Newman, J. N. (1977), Marine Hydrodynamics; Hess and Smith (1964); project boundary-condition specification.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% --- 1. Validate Inputs and Initialize the Algorithm ---

    if nargin < 3 || isempty(factors)
        if size(A, 1) ~= size(A, 2), error('CRESTU:SystemShape', 'A must be square.'); end
        if size(B, 1) ~= size(A, 1), error('CRESTU:RhsShape', 'A and B row counts differ.'); end
        [L, U, p] = lu(A, 'vector');
        factors = struct('L', L, 'U', U, 'p', p, 'n', size(A, 1));
    else
        if size(B, 1) ~= factors.n, error('CRESTU:RhsShape', 'Cached LU and B row counts differ.'); end
        L = factors.L; U = factors.U; p = factors.p;
    end
    X = U \ (L \ B(p, :));
end
