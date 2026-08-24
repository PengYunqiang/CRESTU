function result = solve_radiation_freq(domain,omega,rho,g)
% SOLVE_RADIATION_FREQ Execute the documented solve_radiation_freq operation.
%
% Syntax:
%   result = solve_radiation_freq(domain,omega,rho,g)
%
% Inputs:
%   domain          : [struct] Assembled Rankine boundary domain and configuration.
%   omega           : [scalar] Angular frequency, in rad/s.
%   rho             : [scalar] Water density, in kg/m^3.
%   g               : [scalar] Gravitational acceleration, in m/s^2.
%
% Outputs:
%   result          : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
%SOLVE_RADIATION_FREQ Solve all 6N radiation modes for one frequency.
% This compatibility entry point uses the same assembly/LU/force modules as
% run_frequency_domain_case.
    if nargin<3||isempty(rho), rho=domain.cfg.rho; end
    if nargin<4||isempty(g), g=domain.cfg.grav; end
    if domain.cfg.isx||domain.cfg.isy
        error('CRESTU:SymmetryCompatibilityEntry', ...
            ['solve_radiation_freq is a full-domain compatibility entry; ', ...
             'use run_frequency_domain_case for ISX/ISY cases.']);
    end
    if isfield(domain,'geometry'), geom=domain.geometry; else, geom=merge_domain_geometry(domain); end
    nb=domain.stats.total_body_panels;
    nj=compute_generalized_normals(geom.centers(1:nb,:),geom.normals(1:nb,:),domain.body_list);
    [K,S]=assemble_rankine_matrix(geom.total_panels,geom.centers,geom.normals,geom.vertices, ...
        domain.stats,omega,g,domain.cfg);
    phi=solve_complex_system(K,S*(1i*omega*nj)); phi_body=phi(1:nb,:);
    [A,B,diagnostics]=compute_hydrodynamic_coeffs(phi_body,nj,geom.areas(1:nb),omega,rho);
    result=struct('omega',omega,'phi',phi,'phi_body',phi_body,'added_mass',A, ...
        'damping',B,'generalized_normals',nj,'diagnostics',diagnostics);
end
