function projectRoot = setup_crestu()
%SETUP_CRESTU Add the portable CRESTU source paths.
%   PROJECTROOT = SETUP_CRESTU() resolves paths from this file, not pwd.

    projectRoot = fileparts(mfilename('fullpath'));
    sourceRoot = fullfile(projectRoot, 'Source Code');
    requiredDirectories = {
        fullfile(sourceRoot, '0.Tools Code')
        fullfile(sourceRoot, '1.Input')
        fullfile(sourceRoot, '2.Mesh')
        fullfile(sourceRoot, '3.HessSmith')
        fullfile(sourceRoot, '4.Potential')
        fullfile(sourceRoot, '5.Force')
        fullfile(sourceRoot, '6.MeanDriftLoads')
        fullfile(projectRoot, 'Case_Wave')
        };

    for directoryIndex = 1:numel(requiredDirectories)
        directory = requiredDirectories{directoryIndex};
        assert(isfolder(directory), 'CRESTU:MissingSourceDirectory', ...
            'Required CRESTU directory was not found: %s', directory);
        addpath(directory, '-begin');
    end

    matlabInfo = ver('MATLAB');
    fprintf('[OK] CRESTU paths configured from %s\n', projectRoot);
    fprintf('[INFO] MATLAB %s\n', matlabInfo.Version);
end
