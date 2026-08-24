function [phi_I,dphi_I_dn,k] = compute_incident_wave(centers,normals,omega,g,depth,beta_deg,amplitude)
% COMPUTE_INCIDENT_WAVE Execute the documented compute_incident_wave operation.
%
% Syntax:
%   [phi_I,dphi_I_dn,k] = compute_incident_wave(centers,normals,omega,g,depth,beta_deg,amplitude)
%
% Inputs:
%   centers         : [N x 3] Panel collocation points, in m.
%   normals         : [N x 3] Unit panel normals pointing into the fluid.
%   omega           : [scalar] Angular frequency, in rad/s.
%   g               : [scalar] Gravitational acceleration, in m/s^2.
%   depth           : [scalar] Water depth, in m; zero denotes infinite depth.
%   beta_deg        : [scalar] Incident-wave heading, in deg.
%   amplitude       : [scalar] Incident-wave amplitude, in m.
%
% Outputs:
%   phi_I           : [documented value] Function result; dimensions and units follow the implemented contract.
%   dphi_I_dn       : [documented value] Function result; dimensions and units follow the implemented contract.
%   k               : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
%COMPUTE_INCIDENT_WAVE Airy incident potential for exp(i*omega*t).
% beta is the direction of propagation measured counterclockwise from +x.
    if nargin<7||isempty(amplitude), amplitude=1; end
    centers=reshape(centers,[],3); normals=reshape(normals,[],3);
    if size(centers,1)~=size(normals,1), error('CRESTU:IncidentShape','Centers/normals row mismatch.'); end
    validateattributes(omega,{'numeric'},{'scalar','real','positive','finite'});
    k=solve_dispersion(omega,g,depth); beta=beta_deg*pi/180;
    phase=exp(-1i*k*(centers(:,1)*cos(beta)+centers(:,2)*sin(beta)));
    coefficient=1i*g*amplitude/omega;
    if depth>0
        vertical=cosh(k*(centers(:,3)+depth))/cosh(k*depth);
        vertical_dz=k*sinh(k*(centers(:,3)+depth))/cosh(k*depth);
    else
        vertical=exp(k*centers(:,3)); vertical_dz=k*vertical;
    end
    phi_I=coefficient*vertical.*phase;
    u=(-1i*k*cos(beta))*phi_I; v=(-1i*k*sin(beta))*phi_I;
    w=coefficient*vertical_dz.*phase;
    dphi_I_dn=u.*normals(:,1)+v.*normals(:,2)+w.*normals(:,3);
end

function k=solve_dispersion(omega,g,depth)
% SOLVE_DISPERSION Execute the documented solve_dispersion operation.
%
% Syntax:
%   k=solve_dispersion(omega,g,depth)
%
% Inputs:
%   omega           : [scalar] Angular frequency, in rad/s.
%   g               : [scalar] Gravitational acceleration, in m/s^2.
%   depth           : [scalar] Water depth, in m; zero denotes infinite depth.
%
% Outputs:
%   k               : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
    k=omega^2/g;
    if depth<=0, return; end
    k=max(k,omega/sqrt(g*depth));
    for iter=1:50
        kh=k*depth; t=tanh(kh); f=g*k*t-omega^2;
        df=g*t+g*k*depth/(cosh(kh)^2); step=f/df; k_new=max(k-step,eps);
        if abs(k_new-k)<=1e-12*max(1,k), k=k_new; return; end
        k=k_new;
    end
    warning('CRESTU:DispersionConvergence','Dispersion iteration reached its limit at omega=%g.',omega);
end
