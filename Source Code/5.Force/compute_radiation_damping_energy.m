function [damping,diagnostics] = compute_radiation_damping_energy( ...
        phi_radiation,nj,centers,normals,areas,omega,cfg,mode_parity,n_theta)
% COMPUTE_RADIATION_DAMPING_ENERGY Execute the documented compute_radiation_damping_energy operation.
%
% Syntax:
%   [damping,diagnostics] = compute_radiation_damping_energy( ...
%       phi_radiation,nj,centers,normals,areas,omega,cfg,mode_parity,n_theta)
%
% Inputs:
%   phi_radiation   : [N x Ndof] Complex radiation potential, in m^2/s.
%   nj              : [N x 6Nb] Generalized panel normals for all body modes.
%   centers         : [N x 3] Panel collocation points, in m.
%   normals         : [N x 3] Unit panel normals pointing into the fluid.
%   areas           : [N x 1] Panel areas, in m^2.
%   omega           : [scalar] Angular frequency, in rad/s.
%   cfg             : [struct] Validated CRESTU configuration with SI-valued physical fields.
%   mode_parity     : [documented value] Input required by the implemented function contract.
%   n_theta         : [integer scalar or array] Discrete count or index required by the algorithm.
%
% Outputs:
%   damping         : [Ndof x Ndof x Nf] Radiation-damping matrices in SI translational/rotational units.
%   diagnostics     : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
%COMPUTE_RADIATION_DAMPING_ENERGY Positive damping from Haskind energy flux.
% B_ij=omega*k/(4*pi*rho*g^2*D)*int X_i(beta)conj(X_j(beta))dbeta,
% D=tanh(kh)+kh/cosh(kh)^2.  This avoids extracting a small dissipative
% real part from the otherwise reactive near-field pressure integral.
    if nargin<9||isempty(n_theta), n_theta=72; end
    headings=(0:n_theta-1)*(360/n_theta); dtheta=2*pi/n_theta;
    force=compute_haskind_excitation(phi_radiation,nj,centers,normals,areas, ...
        omega,cfg,headings,mode_parity);
    [wave_number,~]=solve_wave_dispersion(omega,cfg.grav,cfg.water_depth);
    if cfg.water_depth>0
        kh=wave_number*cfg.water_depth;
        depth_factor=tanh(kh)+kh/(cosh(kh)^2);
    else
        depth_factor=1;
    end
    prefactor=omega*wave_number/(4*pi*cfg.rho*cfg.grav^2*depth_factor);
    damping=real(prefactor*(force*force')*dtheta);
    damping=0.5*(damping+damping.');
    eigenvalues=eig(damping); tolerance=1e-10*max(1,max(abs(eigenvalues)));
    if min(eigenvalues)<-tolerance
        warning('CRESTU:EnergyDampingPSD','Energy damping has minimum eigenvalue %g.',min(eigenvalues));
    end
    diagnostics=struct('method','Haskind energy flux','headings',headings, ...
        'haskind_excitation',force,'depth_factor',depth_factor,'wavenumber',wave_number, ...
        'min_eigenvalue',min(eigenvalues));
end
