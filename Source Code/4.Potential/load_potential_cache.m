function pot_cache = load_potential_cache(cache_file,cfg,geometry)
% LOAD_POTENTIAL_CACHE Execute the documented load_potential_cache operation.
%
% Syntax:
%   pot_cache = load_potential_cache(cache_file,cfg,geometry)
%
% Inputs:
%   cache_file      : [char|string] Potential-cache MAT-file path.
%   cfg             : [struct] Validated CRESTU configuration with SI-valued physical fields.
%   geometry        : [struct] Concatenated panel geometry for all boundary components.
%
% Outputs:
%   pot_cache       : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
%LOAD_POTENTIAL_CACHE Load and validate a CRESTU potential cache.
    if ~exist(cache_file,'file'), error('CRESTU:CacheMissing','Potential cache not found: %s',cache_file); end
    data=load(cache_file,'pot_cache');
    if ~isfield(data,'pot_cache'), error('CRESTU:CacheVariable','Cache lacks pot_cache variable: %s',cache_file); end
    pot_cache=data.pot_cache;
    required={'schema_version','case_name','omegas','headings','n_body_panels', ...
        'n_total_panels','environment','geometry_signature','entries'};
    for k=1:numel(required)
        if ~isfield(pot_cache,required{k}), error('CRESTU:CacheSchema','Cache lacks %s.',required{k}); end
    end
    active_environment=[cfg.water_depth,cfg.grav,cfg.rho,cfg.fs.r_inner,cfg.fs.r_outer, ...
        cfg.fs.mu0,cfg.isx,cfg.isy];
    active_signature=geometry_signature(geometry);
    if pot_cache.schema_version~=4||~strcmp(pot_cache.case_name,cfg.case_name) || ...
            pot_cache.n_body_panels~=geometry.body_panels || pot_cache.n_total_panels~=geometry.total_panels || ...
            ~same_array(pot_cache.environment,active_environment) || ...
            ~same_array(pot_cache.geometry_signature,active_signature) || ...
            ~same_array(pot_cache.omegas,cfg.freq.omegas) || ...
            ~same_array(pot_cache.headings,cfg.wave.headings)
        error('CRESTU:CacheMismatch','Cache metadata do not match the active case/configuration.');
    end
end

function signature=geometry_signature(g)
% GEOMETRY_SIGNATURE Execute the documented geometry_signature operation.
%
% Syntax:
%   signature=geometry_signature(g)
%
% Inputs:
%   g               : [scalar] Gravitational acceleration, in m/s^2.
%
% Outputs:
%   signature       : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
    signature=[sum(g.centers(:)),sum(g.centers(:).^2),sum(g.normals(:)), ...
        sum(g.areas(:)),sum(g.vertices(:)),sum(g.vertices(:).^2)];
end

function tf=same_array(a,b)
% SAME_ARRAY Execute the documented same_array operation.
%
% Syntax:
%   tf=same_array(a,b)
%
% Inputs:
%   a               : [documented value] Input required by the implemented function contract.
%   b               : [documented value] Input required by the implemented function contract.
%
% Outputs:
%   tf              : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
    a=a(:); b=b(:); tf=numel(a)==numel(b)&&all(abs(a-b)<=1e-12*max(1,max(abs([a;b]))));
end
