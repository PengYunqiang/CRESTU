function [results,manifest]=CRESTU_2D_BATCH(inputs,outputRoot)
%CRESTU_2D_BATCH Sequential independent sections; no strip-theory coupling.
% [results,manifest]=CRESTU_2D_BATCH({file1,cfg2},newOutputDirectory)
% Each input keeps its frequencies, DOFs and output.save setting. The explicit
% batch root replaces output.directory with case_001, case_002, ... folders.
    CRESTU_SETUP('2d');
    if ischar(inputs) || (isstring(inputs) && isscalar(inputs)) || isstruct(inputs)
        inputs={inputs};
    elseif isstring(inputs)
        inputs=cellstr(inputs);
    end
    assert(iscell(inputs) && isvector(inputs) && ~isempty(inputs),'CRESTU2D:BatchInput', ...
        'Supply a nonempty list of JSON filenames or scalar case structures.');
    assert((ischar(outputRoot) && isrow(outputRoot) && ~isempty(outputRoot)) || ...
        (isstring(outputRoot) && isscalar(outputRoot) && strlength(outputRoot)>0), ...
        'CRESTU2D:Output','Batch outputRoot must be nonempty scalar text.');
    outputRoot=char(outputRoot);
    if ~(startsWith(outputRoot,'/') || startsWith(outputRoot,'\\') || ...
            (numel(outputRoot)>=2 && outputRoot(2)==':'))
        outputRoot=fullfile(pwd,outputRoot);
    end
    assert(~isfolder(outputRoot) && ~isfile(outputRoot),'CRESTU2D:OutputExists', ...
        'Choose a new batch output directory; prior results are retained.');
    cfgs=cell(size(inputs)); results=cell(size(inputs)); entries=struct([]);
    % Preflight every section before starting any hydrodynamic solve.
    for j=1:numel(inputs)
        cfg=crestu2d.read_case(inputs{j}); crestu2d.geometry(cfg);
        cfg.output.directory=fullfile(outputRoot,sprintf('case_%03d',j)); cfgs{j}=cfg;
        entry=struct('index',j,'name',cfg.name,'problem',cfg.problem, ...
            'dofs',{cfg.dofs},'omega_rad_s',cfg.omega_rad_s, ...
            'output_directory',cfg.output.directory,'save_requested',cfg.output.save, ...
            'status','NOT_RUN','elapsed_seconds',[],'error_identifier','','error_message','');
        % Establish the field schema before assigning later indexed entries.
        if j==1,entries=entry;else,entries(j)=entry;end %#ok<AGROW>
    end
    mkdir(outputRoot);
    manifest=struct('schema','CRESTU-2D-BATCH-1','status','RUNNING', ...
        'scope','Sequential independent section BEM solves; no coupling, interpolation or strip integration.', ...
        'started_utc',utc_now(),'output_directory',outputRoot,'entries',entries);
    write_manifest();
    for j=1:numel(cfgs)
        manifest.entries(j).status='RUNNING'; write_manifest();
        try
            results{j}=CRESTU_2D(cfgs{j});
            manifest.entries(j).status='RUNTIME_COMPLETED';
            manifest.entries(j).elapsed_seconds=results{j}.metadata.elapsed_seconds;
            write_manifest();
        catch ex
            manifest.entries(j).status='FAILED';
            manifest.entries(j).error_identifier=ex.identifier;
            manifest.entries(j).error_message=ex.message;
            manifest.status='FAILED'; manifest.completed_utc=utc_now(); write_manifest();
            rethrow(ex);
        end
    end
    manifest.status='COMPLETED_INDEPENDENT_SECTIONS'; manifest.completed_utc=utc_now(); write_manifest();

    function write_manifest()
        fid=fopen(fullfile(outputRoot,'batch_manifest.json'),'w','n','UTF-8');assert(fid>=0);
        cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
        fprintf(fid,'%s\n',jsonencode(manifest,PrettyPrint=true));
    end
end
function text=utc_now()
    text=char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss.SSS''Z'''));
end
