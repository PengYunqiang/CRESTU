function save_potential_cache(cache_file,pot_cache)
% SAVE_POTENTIAL_CACHE Execute the documented save_potential_cache operation.
%
% Syntax:
%   save_potential_cache(cache_file,pot_cache)
%
% Inputs:
%   cache_file      : [char|string] Potential-cache MAT-file path.
%   pot_cache       : [documented value] Input required by the implemented function contract.
%
% Outputs:
%   None; the function performs the documented file, plot, or validation action.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
%SAVE_POTENTIAL_CACHE Save a complete potential sweep with schema metadata.
    if ~isstruct(pot_cache)||~isfield(pot_cache,'schema_version')
        error('CRESTU:CacheSchema','Potential cache must contain schema_version.');
    end
    cache_dir=fileparts(cache_file);
    if ~isempty(cache_dir)&&~exist(cache_dir,'dir'), mkdir(cache_dir); end
    save(cache_file,'pot_cache','-v7.3');
    fprintf('>>> Potential cache saved: %s\n',cache_file);
end
