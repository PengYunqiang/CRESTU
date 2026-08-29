function potentialCache = load_potential_cache(cacheFilename, expectedBaseKey)
% LOAD_POTENTIAL_CACHE Load a potential cache only after exact key validation.
%
% Syntax:
%   potentialCache = load_potential_cache(cacheFilename, expectedBaseKey)
%
% Inputs:
%   cacheFilename   - Potential-cache MAT-file path [-].
%   expectedBaseKey - SHA-256 key for the active case and base geometry [-].
%
% Outputs:
%   potentialCache  - Schema-5 potential cache with per-frequency keys [-].

    arguments
        cacheFilename {mustBeTextScalar}
        expectedBaseKey {mustBeTextScalar}
    end

    %% 阶段 1: 验证缓存文件与变量结构

    if ~isfile(cacheFilename)
        fprintf('[CACHE] MISS | file=%s | reason=file_missing\n', ...
            cacheFilename);
        error('CRESTU:CacheMissing', ...
            'Potential cache not found: %s', cacheFilename);
    end

    loadedData = load(cacheFilename, 'pot_cache');

    if ~isfield(loadedData, 'pot_cache')
        fprintf('[CACHE] MISS | file=%s | reason=variable_missing\n', ...
            cacheFilename);
        error('CRESTU:CacheVariable', ...
            'Cache lacks pot_cache variable: %s', cacheFilename);
    end

    potentialCache = loadedData.pot_cache;
    requiredFields = {'schema_version', 'case_name', 'base_cache_key', ...
        'code_version', 'entries'};

    for fieldIndex = 1:numel(requiredFields)
        if ~isfield(potentialCache, requiredFields{fieldIndex})
            fprintf('[CACHE] MISS | file=%s | reason=schema_field_%s\n', ...
                cacheFilename, requiredFields{fieldIndex});
            error('CRESTU:CacheSchema', 'Cache lacks %s.', ...
                requiredFields{fieldIndex});
        end
    end

    %% 阶段 2: 对比完整基础缓存键

    if potentialCache.schema_version ~= 5
        fprintf('[CACHE] MISS | file=%s | reason=schema_%d\n', ...
            cacheFilename, potentialCache.schema_version);
        error('CRESTU:CacheSchema', ...
            'Potential cache schema %d is obsolete; schema 5 is required.', ...
            potentialCache.schema_version);
    end

    if ~strcmp(potentialCache.base_cache_key, char(expectedBaseKey))
        fprintf('[CACHE] MISS | file=%s | reason=base_key_mismatch\n', ...
            cacheFilename);
        error('CRESTU:CacheMismatch', ...
            'Potential cache base key does not match the active solver state.');
    end

    fprintf('[CACHE] HIT | file=%s | base_key=%s\n', ...
        cacheFilename, potentialCache.base_cache_key);
end
