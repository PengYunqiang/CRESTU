function [indices, matched_omegas, error_value] = match_reference_grid(reference_omegas, target_omegas, tolerance)
% MATCH_REFERENCE_GRID Match requested frequencies to a WAMIT frequency grid.
%
% Syntax:
%   [indices, matched_omegas, error_value] = match_reference_grid(reference_omegas, target_omegas, tolerance)
%
% Inputs:
%   reference_omegas : [1 x Nr] Available angular frequencies, in rad/s.
%   target_omegas    : [1 x Nt] Requested angular frequencies, in rad/s.
%   tolerance        : [scalar] Maximum absolute matching error, in rad/s.
%
% Outputs:
%   indices          : [1 x Nt] Indices into the reference grid.
%   matched_omegas   : [1 x Nt] Actual matched reference frequencies, in rad/s.
%   error_value      : [1 x Nt] Absolute matching errors, in rad/s.
%
% Mathematical Reference:
%   Nearest-neighbor matching on a one-dimensional monotone frequency grid.
    if nargin < 3 || isempty(tolerance)
        tolerance = 2.0e-4;
    end
    reference_omegas = reshape(reference_omegas, 1, []);
    target_omegas = reshape(target_omegas, 1, []);
    distance = abs(reference_omegas.' - target_omegas);
    [error_value, indices] = min(distance, [], 1);
    if any(error_value > tolerance)
        error('CRESTU:ReferenceGrid', 'Reference frequency mismatch exceeds %.3g rad/s.', tolerance);
    end
    matched_omegas = reference_omegas(indices);
end
