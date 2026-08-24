function force = compute_haskind_excitation(phi_radiation,nj,centers,normals,areas,omega,cfg,headings,mode_parity)
% COMPUTE_HASKIND_EXCITATION Execute the documented compute_haskind_excitation operation.
%
% Syntax:
%   force = compute_haskind_excitation(phi_radiation,nj,centers,normals,areas,omega,cfg,headings,mode_parity)
%
% Inputs:
%   phi_radiation   : [N x Ndof] Complex radiation potential, in m^2/s.
%   nj              : [N x 6Nb] Generalized panel normals for all body modes.
%   centers         : [N x 3] Panel collocation points, in m.
%   normals         : [N x 3] Unit panel normals pointing into the fluid.
%   areas           : [N x 1] Panel areas, in m^2.
%   omega           : [scalar] Angular frequency, in rad/s.
%   cfg             : [struct] Validated CRESTU configuration with SI-valued physical fields.
%   headings        : [1 x Nh] Incident-wave headings, in deg.
%   mode_parity     : [documented value] Input required by the implemented function contract.
%
% Outputs:
%   force           : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
%COMPUTE_HASKIND_EXCITATION Excitation from radiation potentials.
% Uses unit-velocity radiation potential psi=phi_R/(i*omega):
% With CRESTU normals directed from body into fluid,
% X_j=-i*omega*rho*int[n_j*phi_I+psi_j*d(phi_I)/dn]dS.
    areas=reshape(areas,[],1); headings=reshape(headings,1,[]);
    ndof=size(nj,2); nh=numel(headings); force=complex(zeros(ndof,nh));
    if nargin<9||isempty(mode_parity)
        mode_parity=get_mode_parities(cfg.n_bodies,cfg.isx,cfg.isy);
    end
    psi=phi_radiation/(1i*omega);
    for j=1:ndof
        parity=mode_parity(j,:);
        for h=1:nh
            [phi_I,dphi_I]=decompose_incident_wave_symmetry(centers,normals,omega,cfg.grav, ...
                cfg.water_depth,headings(h),1,cfg.isx,cfg.isy,parity);
            integral=sum((nj(:,j).*phi_I+psi(:,j).*dphi_I).*areas);
            force(j,h)=-1i*omega*cfg.rho*cfg.symmetry.multiplicity*integral;
        end
    end
end
