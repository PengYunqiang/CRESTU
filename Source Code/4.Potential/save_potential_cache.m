function save_potential_cache(cache_file, pot_cache)
% SAVE_POTENTIAL_CACHE Save potential cache for the CRESTU hydrodynamic workflow.
%
% Syntax:
%   save_potential_cache(cache_file, pot_cache)
%
% Description:
%   Implements linear Rankine potential-flow operations.
%   Symmetry and phase follow exp(i*omega*t).
%
% Inputs:
%   cache_file         - [character vector or string scalar] Potential-cache file path.
%   pot_cache          - [struct] Cached geometry signature, frequency metadata, and complex potentials.
%
% Outputs:
%   None.
%
% Governing Equations / Theory:
%   Green third identity, the Rankine kernel 1/r, linearized free-surface theory, and reflection symmetry as applicable.
%
% References:
%   - Newman, J. N. (1977), Marine Hydrodynamics; Hess and Smith (1964); project boundary-condition specification.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    if ~isstruct(pot_cache) || ~isfield(pot_cache,'schema_version')
        error('CRESTU:CacheSchema','Potential cache must contain schema_version.');
    end
    cache_dir = fileparts(cache_file);
    if ~isempty(cache_dir) && ~exist(cache_dir,'dir')
        mkdir(cache_dir);
    end
    save(cache_file,'pot_cache','-v7.3');
    fprintf('[OK] Potential cache saved: %s\n', cache_file);
end
