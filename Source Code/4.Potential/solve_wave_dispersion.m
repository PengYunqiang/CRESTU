function [k, wavelength] = solve_wave_dispersion(omega, g, depth)
% SOLVE_WAVE_DISPERSION Solve the linear gravity-wave dispersion relation for the wavenumber.
%
% Syntax:
%   [k, wavelength] = solve_wave_dispersion(omega, g, depth)
%
% Description:
%   Implements linear Rankine potential-flow operations.
%   Symmetry and phase follow exp(i*omega*t).
%
% Inputs:
%   omega              - [scalar] Angular frequency, [rad/s].
%   g                  - [scalar] Gravitational acceleration, [m/s^2].
%   depth              - [scalar] Water depth; a nonpositive value denotes infinite depth, [m].
%
% Outputs:
%   k                  - [scalar] Gravity-wave wavenumber, [1/m].
%   wavelength         - [scalar] Gravity wavelength, [m].
%
% Governing Equations / Theory:
%   Green third identity, the Rankine kernel 1/r, linearized free-surface theory, and reflection symmetry as applicable.
%
% References:
%   - Newman, J. N. (1977), Marine Hydrodynamics; Hess and Smith (1964); project boundary-condition specification.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    k = omega^2 / g;
    if depth > 0
        k = max(k, omega / sqrt(g * depth));
        for iter = 1:50
            kh = k * depth;
            t = tanh(kh);
            dispersionResidual = g * k * t - omega^2;
            df = g * t + g * k * depth / (cosh(kh)^2);
            candidate = max(k - dispersionResidual / df, eps);
            if abs(candidate - k) <= 1e-12 * max(1, k)
                k = candidate;
                break;
            end
            k = candidate;
        end
    end
    wavelength = 2 * pi / k;
end
