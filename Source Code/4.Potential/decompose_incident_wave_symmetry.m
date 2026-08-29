function [phi_component, dphi_component] = decompose_incident_wave_symmetry( ...
    centers, normals, omega, g, depth, beta_deg, amplitude, isx, isy, parity)
% DECOMPOSE_INCIDENT_WAVE_SYMMETRY Decompose incident-wave boundary data into reflection-parity components.
%
% Syntax:
%   [phi_component, dphi_component] = decompose_incident_wave_symmetry(centers, normals, omega, g, depth, beta_deg, amplitude, isx, isy, parity)
%
% Description:
%   Implements linear Rankine potential-flow operations.
%   Symmetry and phase follow exp(i*omega*t).
%
% Inputs:
%   centers            - [N x 3] Panel collocation points in global coordinates, [m].
%   normals            - [N x 3] Unit panel normals, dimensionless.
%   omega              - [scalar] Angular frequency, [rad/s].
%   g                  - [scalar] Gravitational acceleration, [m/s^2].
%   depth              - [scalar] Water depth; a nonpositive value denotes infinite depth, [m].
%   beta_deg           - [scalar or vector] Wave-propagation heading, [deg].
%   amplitude          - [scalar] Incident-wave amplitude, [m].
%   isx                - [logical scalar] Symmetry-reduction flag for the x = 0 plane.
%   isy                - [logical scalar] Symmetry-reduction flag for the y = 0 plane.
%   parity             - [1 x 2] Reflection parity signs for the x = 0 and y = 0 planes, dimensionless.
%
% Outputs:
%   phi_component      - [N x Nh] Incident-potential component with the requested symmetry parity, [m^2/s].
%   dphi_component     - [N x Nh] Normal derivative of the parity-decomposed incident potential, [m/s].
%
% Governing Equations / Theory:
%   Green third identity, the Rankine kernel 1/r, linearized free-surface theory, and reflection symmetry as applicable.
%
% References:
%   - Newman, J. N. (1977), Marine Hydrodynamics; Hess and Smith (1964); project boundary-condition specification.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    flags = [0, 0];
    headings = beta_deg;
    weights = 1;
    if isx
        flags = [flags;1, 0];
        headings = [headings, 180 - beta_deg];
        weights = [weights, parity(1)];
    end
    if isy
        flags = [flags;0, 1];
        headings = [headings, -beta_deg];
        weights = [weights, parity(2)];
    end
    if isx && isy
        flags = [0, 0;1, 0;0, 1;1, 1];
        headings = [beta_deg, 180 - beta_deg, -beta_deg, 180 + beta_deg];
        weights = [1, parity(1), parity(2), parity(1) * parity(2)];
    end
    phi_component = complex(zeros(size(centers, 1), 1));
    dphi_component = complex(zeros(size(centers, 1), 1));
    for componentIndex = 1:size(flags, 1)
        [phi, dphi] = compute_incident_wave(centers, normals, omega, g, depth, headings(componentIndex), amplitude);
        phi_component = phi_component + weights(componentIndex) * phi;
        dphi_component = dphi_component + weights(componentIndex) * dphi;
    end
    scale = 2^(isx + isy);
    phi_component = phi_component / scale;
    dphi_component = dphi_component / scale;
end
