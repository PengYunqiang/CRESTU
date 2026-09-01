function remove_runtime(definition, runtimeDirectory)
% REMOVE_RUNTIME Remove only a verified Phase 3.2 isolated runtime folder.

    if ~isfolder(runtimeDirectory)
        return
    end
    artifactRoot = char(java.io.File(definition.artifactDirectory).getCanonicalPath());
    runtimeRoot = char(java.io.File(runtimeDirectory).getCanonicalPath());
    expectedPrefix = [artifactRoot, filesep];
    [~, runtimeName] = fileparts(runtimeRoot);
    assert(startsWith(lower(runtimeRoot), lower(expectedPrefix)) && ...
        startsWith(runtimeName, 'tp'), ...
        'CRESTU:Phase32RuntimeCleanupScope', ...
        'Refusing to remove a directory outside the isolated runtime scope.');
    [removed, message] = rmdir(runtimeRoot, 's');
    assert(removed, 'CRESTU:Phase32RuntimeCleanup', ...
        'Cannot remove isolated runtime directory: %s', message);
end
