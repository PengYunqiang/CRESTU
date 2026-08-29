function [solution, factors] = solve_complex_system(A, B, factors)
% SOLVE_COMPLEX_SYSTEM Solve a reusable complex linear system with cached factorization.
%
% Syntax:
%   [solution, factors] = solve_complex_system(A, B, factors)
%
% Description:
%   Implements linear Rankine potential-flow operations.
%   Symmetry and phase follow exp(i*omega*t).
%
% Inputs:
%   A                  - [numeric array] Matrix or geometric area defined by the governing formulation, in corresponding SI units.
%   B                  - [numeric array] Secondary matrix or geometric breadth defined by the function contract, in corresponding SI units.
%   factors            - [struct or empty] Reusable complex-system factorization data.
%
% Outputs:
%   solution           - [numeric array] Complex solution matrix with units inherited from the right-hand side.
%   factors            - [struct or empty] Reusable complex-system factorization data.
%
% Governing Equations / Theory:
%   Green third identity, the Rankine kernel 1/r, linearized free-surface theory, and reflection symmetry as applicable.
%
% References:
%   - Newman, J. N. (1977), Marine Hydrodynamics; Hess and Smith (1964); project boundary-condition specification.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    % <<<CORE>>> solve_rankine_complex_system, paper_eq=none, benchmark=single_sphere_linear_residual
    if nargin < 3 || isempty(factors)
        if size(A, 1) ~= size(A, 2)
            error('CRESTU:SystemShape','A must be square.');
        end
        if size(B, 1) ~= size(A, 1)
            error('CRESTU:RhsShape','A and B row counts differ.');
        end
        [lowerFactor, upperFactor, permutationVector] = lu(A, 'vector');
        factors = struct('L', lowerFactor, 'U', upperFactor, ...
            'p', permutationVector, 'n', size(A, 1));
    else
        if size(B, 1) ~= factors.n
            error('CRESTU:RhsShape','Cached LU and B row counts differ.');
        end
        lowerFactor = factors.L;
        upperFactor = factors.U;
        permutationVector = factors.p;
    end
    solution = upperFactor \ (lowerFactor \ B(permutationVector, :));
    % <<</CORE>>>
end
