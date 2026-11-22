function root=CRESTU_SETUP(mode)
%CRESTU_SETUP Add the standalone 2-D package and verify source resolution.
    if nargin<1,mode='2d';end
    assert(strcmpi(mode,'2d'),'CRESTU:Dimension','This package supports mode 2d.');
    root=fileparts(mfilename('fullpath'));addpath(root);
    names={'CRESTU_SETUP','CRESTU_2D','CRESTU_2D_BATCH','CRESTU_2D_HF'};
    expected=cellfun(@(n)fullfile(root,[n '.m']),names,'UniformOutput',false);
    packageRoot=fullfile(root,'source','2d');addpath(packageRoot);
    files=dir(fullfile(packageRoot,'+crestu2d','*.m'));
    for j=1:numel(files)
        [~,name]=fileparts(files(j).name);
        names{end+1}=['crestu2d.' name]; %#ok<AGROW>
        expected{end+1}=fullfile(files(j).folder,files(j).name); %#ok<AGROW>
    end
    for j=1:numel(names)
        resolved=which(names{j});[okA,a]=fileattrib(expected{j});[okB,b]=fileattrib(resolved);
        same=false;
        if okA && okB
            if ispc,same=strcmpi(a.Name,b.Name);else,same=strcmp(a.Name,b.Name);end
        end
        assert(same,'CRESTU:SourcePathShadowing', ...
            'Function %s must resolve to %s; actual %s. Use one CRESTU tree on the MATLAB path.', ...
            names{j},expected{j},resolved);
    end
end
