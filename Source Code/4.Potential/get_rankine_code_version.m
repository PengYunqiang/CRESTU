function codeVersion = get_rankine_code_version()
% GET_RANKINE_CODE_VERSION Fingerprint ALL numerical MATLAB source dependencies.
% Includes geometry, incident waves, generalized normals, Haskind, hydrostatics,
% mass, readers and this function. Conservative invalidation is intentional.
    potentialDirectory = fileparts(mfilename('fullpath'));
    sourceDirectory = fileparts(potentialDirectory);
    projectDirectory = fileparts(sourceDirectory);
    listing = dir(fullfile(sourceDirectory, '**', '*.m'));
    relativeFiles = cell(1,numel(listing)+1);
    for k = 1:numel(listing)
        absolute = fullfile(listing(k).folder,listing(k).name);
        relativeFiles{k} = absolute(numel(projectDirectory)+2:end);
    end
    relativeFiles{end} = fullfile('Case_Wave','run_frequency_domain_case.m');
    relativeFiles = sort(relativeFiles);
    fileHashes = cell(size(relativeFiles));
    for k = 1:numel(relativeFiles)
        fileHashes{k} = sha256_file_hash(fullfile(projectDirectory,relativeFiles{k}));
    end
    semanticVersion = 'rankine-gpt6pro-planar-haskind-resolution-v2';
    combinedText = strjoin([{semanticVersion},relativeFiles,fileHashes],'|');
    codeVersion = struct('semanticVersion',semanticVersion, ...
        'fingerprint',sha256_hash(combinedText),'relativeFiles',{relativeFiles}, ...
        'fileHashes',{fileHashes});
end
