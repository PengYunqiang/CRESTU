function [X,factors] = solve_complex_system(A,B,factors)
% SOLVE_COMPLEX_SYSTEM Execute the documented solve_complex_system operation.
%
% Syntax:
%   [X,factors] = solve_complex_system(A,B,factors)
%
% Inputs:
%   A               : [documented value] Input required by the implemented function contract.
%   B               : [documented value] Input required by the implemented function contract.
%   factors         : [documented value] Input required by the implemented function contract.
%
% Outputs:
%   X               : [documented value] Function result; dimensions and units follow the implemented contract.
%   factors         : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
%SOLVE_COMPLEX_SYSTEM LU factorization/reuse for dense complex systems.
    if nargin<3||isempty(factors)
        if size(A,1)~=size(A,2), error('CRESTU:SystemShape','A must be square.'); end
        if size(B,1)~=size(A,1), error('CRESTU:RhsShape','A and B row counts differ.'); end
        [L,U,p]=lu(A,'vector');
        factors=struct('L',L,'U',U,'p',p,'n',size(A,1));
    else
        if size(B,1)~=factors.n, error('CRESTU:RhsShape','Cached LU and B row counts differ.'); end
        L=factors.L; U=factors.U; p=factors.p;
    end
    X=U\(L\B(p,:));
end
