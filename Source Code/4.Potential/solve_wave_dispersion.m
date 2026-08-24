function [k,wavelength] = solve_wave_dispersion(omega,g,depth)
% SOLVE_WAVE_DISPERSION Execute the documented solve_wave_dispersion operation.
%
% Syntax:
%   [k,wavelength] = solve_wave_dispersion(omega,g,depth)
%
% Inputs:
%   omega           : [scalar] Angular frequency, in rad/s.
%   g               : [scalar] Gravitational acceleration, in m/s^2.
%   depth           : [scalar] Water depth, in m; zero denotes infinite depth.
%
% Outputs:
%   k               : [documented value] Function result; dimensions and units follow the implemented contract.
%   wavelength      : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
%SOLVE_WAVE_DISPERSION Solve omega^2=g*k*tanh(k*h); depth<=0 is deep water.
    k=omega^2/g;
    if depth>0
        k=max(k,omega/sqrt(g*depth));
        for iter=1:50
            kh=k*depth; t=tanh(kh); f=g*k*t-omega^2;
            df=g*t+g*k*depth/(cosh(kh)^2); candidate=max(k-f/df,eps);
            if abs(candidate-k)<=1e-12*max(1,k), k=candidate; break; end
            k=candidate;
        end
    end
    wavelength=2*pi/k;
end
