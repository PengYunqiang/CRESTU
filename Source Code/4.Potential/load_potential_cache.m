function pot_cache = load_potential_cache(cache_file, cfg, geometry)
% LOAD_POTENTIAL_CACHE Load potential cache for the CRESTU hydrodynamic workflow.
%
% Syntax:
%   pot_cache = load_potential_cache(cache_file, cfg, geometry)
%
% Description:
%   The routine implements a component of the linear Rankine boundary-element formulation for incompressible, irrotational gravity-wave flow. Geometry, reflection parity, free-surface impedance, and complex phase follow the project convention exp(i*omega*t).
%
% Inputs:
%   cache_file         - [character vector or string scalar] Potential-cache file path.
%   cfg                - [struct] Validated CRESTU configuration containing SI-valued physical and numerical parameters.
%   geometry           - [struct] Merged boundary geometry with coordinates in [m] and areas in [m^2].
%
% Outputs:
%   pot_cache          - [struct] Validated cached potential solution and geometry signature.
%
% Governing Equations / Theory:
%   Green third identity, the Rankine kernel 1/r, linearized free-surface theory, and reflection symmetry as applicable.
%
% References:
%   - Newman, J. N. (1977), Marine Hydrodynamics; Hess and Smith (1964); project boundary-condition specification.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% --- 1. Validate Inputs and Initialize the Algorithm ---

    if ~exist(cache_file, 'file'), error('CRESTU:CacheMissing', 'Potential cache not found: %s', cache_file); end
    data = load(cache_file, 'pot_cache');
    if ~isfield(data, 'pot_cache'), error('CRESTU:CacheVariable', 'Cache lacks pot_cache variable: %s', cache_file); end
    pot_cache = data.pot_cache;
    required = {'schema_version', 'case_name', 'omegas', 'headings', 'n_body_panels', ...
        'n_total_panels', 'environment', 'geometry_signature', 'entries'};
    for k = 1:numel(required)
        if ~isfield(pot_cache, required{k}), error('CRESTU:CacheSchema', 'Cache lacks %s.', required{k}); end
    end
    active_environment = [cfg.water_depth, cfg.grav, cfg.rho, cfg.fs.r_inner, cfg.fs.r_outer, ...
        cfg.fs.mu0, cfg.isx, cfg.isy];
    active_signature = geometry_signature(geometry);
    if pot_cache.schema_version ~= 4 || ~strcmp(pot_cache.case_name, cfg.case_name) || ...
            pot_cache.n_body_panels ~= geometry.body_panels || pot_cache.n_total_panels ~= geometry.total_panels || ...
            ~same_array(pot_cache.environment, active_environment) || ...
            ~same_array(pot_cache.geometry_signature, active_signature) || ...
            ~same_array(pot_cache.omegas, cfg.freq.omegas) || ...
            ~same_array(pot_cache.headings, cfg.wave.headings)
        error('CRESTU:CacheMismatch', 'Cache metadata do not match the active case/configuration.');
    end
end

function signature = geometry_signature(g)
% GEOMETRY_SIGNATURE Construct a deterministic geometry signature for cache validation.
%
% Syntax:
%   signature = geometry_signature(g)
%
% Description:
%   The routine implements a component of the linear Rankine boundary-element formulation for incompressible, irrotational gravity-wave flow. Geometry, reflection parity, free-surface impedance, and complex phase follow the project convention exp(i*omega*t).
%
% Inputs:
%   g                  - [scalar] Gravitational acceleration, [m/s^2].
%
% Outputs:
%   signature          - [struct] Deterministic geometry counts, coordinates, and area signature in SI units.
%
% Governing Equations / Theory:
%   Green third identity, the Rankine kernel 1/r, linearized free-surface theory, and reflection symmetry as applicable.
%
% References:
%   - Newman, J. N. (1977), Marine Hydrodynamics; Hess and Smith (1964); project boundary-condition specification.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% --- 1. Validate Inputs and Initialize the Algorithm ---

    signature = [sum(g.centers(:)), sum(g.centers(:) .^ 2), sum(g.normals(:)), ...
        sum(g.areas(:)), sum(g.vertices(:)), sum(g.vertices(:) .^ 2)];
end

function tf = same_array(a, b)
% SAME_ARRAY Compare two numeric arrays using a scale-aware tolerance.
%
% Syntax:
%   tf = same_array(a, b)
%
% Description:
%   The routine implements a component of the linear Rankine boundary-element formulation for incompressible, irrotational gravity-wave flow. Geometry, reflection parity, free-surface impedance, and complex phase follow the project convention exp(i*omega*t).
%
% Inputs:
%   a                  - First vector or scalar operand; dimensions and SI units follow the stated algorithm.
%   b                  - Second vector or scalar operand; dimensions and SI units follow the stated algorithm.
%
% Outputs:
%   tf                 - [logical scalar] True when the tested condition is satisfied.
%
% Governing Equations / Theory:
%   Green third identity, the Rankine kernel 1/r, linearized free-surface theory, and reflection symmetry as applicable.
%
% References:
%   - Newman, J. N. (1977), Marine Hydrodynamics; Hess and Smith (1964); project boundary-condition specification.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% --- 1. Validate Inputs and Initialize the Algorithm ---

    a = a(:); b = b(:); tf = numel(a) == numel(b) && all(abs(a - b) <= 1e-12 * max(1, max(abs([a;b]))));
end
