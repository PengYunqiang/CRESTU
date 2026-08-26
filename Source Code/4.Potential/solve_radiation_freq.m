function result = solve_radiation_freq(domain, omega, rho, g)
% SOLVE_RADIATION_FREQ Solve all rigid-body radiation modes at one angular frequency.
%
% Syntax:
%   result = solve_radiation_freq(domain, omega, rho, g)
%
% Description:
%   The routine implements a component of the linear Rankine boundary-element formulation for incompressible, irrotational gravity-wave flow. Geometry, reflection parity, free-surface impedance, and complex phase follow the project convention exp(i*omega*t).
%
% Inputs:
%   domain             - [struct] Assembled body, free-surface, seabed, and far-field boundary domain in SI units.
%   omega              - [scalar] Angular frequency, [rad/s].
%   rho                - [scalar] Fluid density, [kg/m^3].
%   g                  - [scalar] Gravitational acceleration, [m/s^2].
%
% Outputs:
%   result             - [struct] Frequency-domain radiation solution and hydrodynamic coefficients in SI units.
%
% Governing Equations / Theory:
%   Green third identity, the Rankine kernel 1/r, linearized free-surface theory, and reflection symmetry as applicable.
%
% References:
%   - Newman, J. N. (1977), Marine Hydrodynamics; Hess and Smith (1964); project boundary-condition specification.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% --- 1. Validate Inputs and Initialize the Algorithm ---

    if nargin < 3 || isempty(rho), rho = domain.cfg.rho; end
    if nargin < 4 || isempty(g), g = domain.cfg.grav; end
    if domain.cfg.isx || domain.cfg.isy
        error('CRESTU:SymmetryCompatibilityEntry', ...
            ['solve_radiation_freq is a full-domain compatibility entry; ', ...
             'use run_frequency_domain_case for ISX/ISY cases.']);
    end
    if isfield(domain, 'geometry'), geom = domain.geometry; else, geom = merge_domain_geometry(domain); end
    nb = domain.stats.total_body_panels;
    nj = compute_generalized_normals(geom.centers(1:nb, :), geom.normals(1:nb, :), domain.body_list);
    [K, S] = assemble_rankine_matrix(geom.total_panels, geom.centers, geom.normals, geom.vertices, ...
        domain.stats, omega, g, domain.cfg);
    phi = solve_complex_system(K, S * (1i * omega * nj)); phi_body = phi(1:nb, :);
    [A, B, diagnostics] = compute_hydrodynamic_coeffs(phi_body, nj, geom.areas(1:nb), omega, rho);
    result = struct('omega', omega, 'phi', phi, 'phi_body', phi_body, 'added_mass', A, ...
        'damping', B, 'generalized_normals', nj, 'diagnostics', diagnostics);
end
