function [phi_component,dphi_component] = decompose_incident_wave_symmetry( ...
    centers,normals,omega,g,depth,beta_deg,amplitude,isx,isy,parity)
% DECOMPOSE_INCIDENT_WAVE_SYMMETRY Execute the documented decompose_incident_wave_symmetry operation.
%
% Syntax:
%   [phi_component,dphi_component] = decompose_incident_wave_symmetry( ...
%       centers,normals,omega,g,depth,beta_deg,amplitude,isx,isy,parity)
%
% Inputs:
%   centers         : [N x 3] Panel collocation points, in m.
%   normals         : [N x 3] Unit panel normals pointing into the fluid.
%   omega           : [scalar] Angular frequency, in rad/s.
%   g               : [scalar] Gravitational acceleration, in m/s^2.
%   depth           : [scalar] Water depth, in m; zero denotes infinite depth.
%   beta_deg        : [scalar] Incident-wave heading, in deg.
%   amplitude       : [scalar] Incident-wave amplitude, in m.
%   isx             : [logical scalar] Reflection-symmetry flag for the x = 0 plane.
%   isy             : [logical scalar] Reflection-symmetry flag for the y = 0 plane.
%   parity          : [K x 2] Reflection parity signs for x and y symmetry planes.
%
% Outputs:
%   phi_component   : [documented value] Function result; dimensions and units follow the implemented contract.
%   dphi_component  : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
%DECOMPOSE_INCIDENT_WAVE_SYMMETRY Project an incident wave onto one parity.
    flags=[0,0]; headings=beta_deg; weights=1;
    if isx
        flags=[flags;1,0]; headings=[headings,180-beta_deg]; weights=[weights,parity(1)];
    end
    if isy
        flags=[flags;0,1]; headings=[headings,-beta_deg]; weights=[weights,parity(2)];
    end
    if isx&&isy
        flags=[0,0;1,0;0,1;1,1];
        headings=[beta_deg,180-beta_deg,-beta_deg,180+beta_deg];
        weights=[1,parity(1),parity(2),parity(1)*parity(2)];
    end
    phi_component=complex(zeros(size(centers,1),1));
    dphi_component=complex(zeros(size(centers,1),1));
    for q=1:size(flags,1)
        [phi,dphi]=compute_incident_wave(centers,normals,omega,g,depth,headings(q),amplitude);
        phi_component=phi_component+weights(q)*phi;
        dphi_component=dphi_component+weights(q)*dphi;
    end
    scale=2^(isx+isy); phi_component=phi_component/scale; dphi_component=dphi_component/scale;
end
