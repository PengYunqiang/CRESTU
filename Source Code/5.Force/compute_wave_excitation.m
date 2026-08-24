function [force,pressure] = compute_wave_excitation(phi_incident,phi_diffraction,nj,areas,omega,rho)
% COMPUTE_WAVE_EXCITATION Execute the documented compute_wave_excitation operation.
%
% Syntax:
%   [force,pressure] = compute_wave_excitation(phi_incident,phi_diffraction,nj,areas,omega,rho)
%
% Inputs:
%   phi_incident    : [N x Nh] Complex incident-wave potential, in m^2/s.
%   phi_diffraction : [N x Nh] Complex diffraction potential, in m^2/s.
%   nj              : [N x 6Nb] Generalized panel normals for all body modes.
%   areas           : [N x 1] Panel areas, in m^2.
%   omega           : [scalar] Angular frequency, in rad/s.
%   rho             : [scalar] Water density, in kg/m^3.
%
% Outputs:
%   force           : [documented value] Function result; dimensions and units follow the implemented contract.
%   pressure        : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
%COMPUTE_WAVE_EXCITATION Integrate FK plus diffraction pressure/force.
% Uses exp(i*omega*t): p=-i*omega*rho*(phi_I+phi_D).
    validateattributes(omega,{'numeric'},{'scalar','real','positive','finite'});
    validateattributes(rho,{'numeric'},{'scalar','real','positive','finite'});
    areas=reshape(areas,[],1);
    if isempty(phi_incident), phi_incident=zeros(size(phi_diffraction),'like',phi_diffraction); end
    if isempty(phi_diffraction), phi_diffraction=zeros(size(phi_incident),'like',phi_incident); end
    if ~isequal(size(phi_incident),size(phi_diffraction)) || ...
            size(phi_incident,1)~=size(nj,1) || numel(areas)~=size(nj,1)
        error('CRESTU:ExcitationShape','Incident/diffraction potentials, normals, and areas are inconsistent.');
    end
    pressure=-1i*omega*rho*(phi_incident+phi_diffraction);
    force=nj.'*(pressure.*areas);
end
