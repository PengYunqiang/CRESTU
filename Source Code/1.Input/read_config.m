function cfg = read_config(config_file)
% READ_CONFIG Execute the documented read_config operation.
%
% Syntax:
%   cfg = read_config(config_file)
%
% Inputs:
%   config_file     : [char|string] Configuration-file path.
%
% Outputs:
%   cfg             : [struct] Validated CRESTU configuration with SI-valued physical fields.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
%READ_CONFIG Parse and validate a CRESTU frequency-domain configuration.
% Common Unicode minus/dash/space characters are accepted. Relative file
% names are resolved against the configuration directory.
    if nargin < 1 || isempty(config_file), config_file = 'CRESTU.cfg'; end
    config_file = absolute_path(config_file);
    fid = fopen(config_file, 'rt');
    if fid < 0, error('CRESTU:ConfigOpen','Cannot open configuration file: %s',config_file); end
    cleanup = onCleanup(@() fclose(fid));
    lines = cell(0,1);
    while ~feof(fid)
        raw = fgetl(fid);
        if ~ischar(raw), continue; end
        raw = strtrim(normalize_text(raw));
        if isempty(raw) || startsWith(raw,'#') || startsWith(raw,'%'), continue; end
        lines{end+1,1} = raw; %#ok<AGROW>
    end
    clear cleanup
    sections = split_sections(lines);
    required = {'PARA1','PARA2','PARA3','PARA4','PARA5','PARA6','PARA7','PARA8','PARA9','PARA10'};
    for k=1:numel(required)
        if ~isfield(sections,required{k})
            error('CRESTU:MissingConfigSection','Missing %s in %s.',required{k},config_file);
        end
    end

    cfg = struct('config_file',config_file,'config_dir',fileparts(config_file));
    cfg.case_name = strtrim(require_line(sections.PARA1,1,'PARA1'));
    if isempty(cfg.case_name), error('CRESTU:InvalidCaseName','PARA1 case name is empty.'); end

    v = numeric_row(require_line(sections.PARA2,1,'PARA2'),'PARA2');
    if numel(v)~=4&&numel(v)~=5
        error('CRESTU:ConfigArrayShape','PARA2 requires 4 flags, or 5 including IDRIFT.');
    end
    validate_flags(v,'PARA2');
    idrift=0;
    if numel(v)==5, idrift=v(5); end
    cfg.run = struct('ipoten',v(1),'iforce',v(2),'irad',v(3),'idiff',v(4), ...
        'idrift',idrift);

    v = numeric_row(require_line(sections.PARA3,1,'PARA3'),'PARA3 frequency line');
    if numel(v)~=4 && numel(v)~=5
        error('CRESTU:ConfigArrayShape','PARA3 frequency line requires 4 or 5 values.');
    end
    cfg.freq.type=v(1); cfg.freq.nfreq=v(2);
    if cfg.freq.nfreq<1 || fix(cfg.freq.nfreq)~=cfg.freq.nfreq
        error('CRESTU:InvalidFrequencyCount','NFREQ must be a positive integer.');
    end
    base=linspace(v(3),v(4),cfg.freq.nfreq);
    expected_step=iff(cfg.freq.nfreq>1,(v(4)-v(3))/(cfg.freq.nfreq-1),0);
    cfg.freq.step=iff(numel(v)==5,v(end),expected_step);
    if numel(v)==5 && cfg.freq.nfreq>1 && abs(cfg.freq.step-expected_step)>1e-10*max(1,abs(expected_step))
        warning('CRESTU:FrequencyStepMismatch', ...
            'FREQ_STEP (%g) conflicts with NFREQ/start/end (%g); NFREQ/start/end are authoritative.', ...
            cfg.freq.step,expected_step);
    end
    switch cfg.freq.type
        case 1, cfg.freq.omegas=base;
        case 2
            if any(base<=0), error('CRESTU:InvalidPeriods','Periods must be positive.'); end
            cfg.freq.omegas=2*pi./base;
        case 3, cfg.freq.omegas=2*pi.*base;
        otherwise, error('CRESTU:InvalidFrequencyType','FREQ_TYPE must be 1, 2, or 3.');
    end
    if any(~isfinite(cfg.freq.omegas)|cfg.freq.omegas<=0)
        error('CRESTU:InvalidFrequencies','All angular frequencies must be finite and positive.');
    end
    cfg.wave.ndir=scalar_integer(require_line(sections.PARA3,2,'PARA3'),'NDIR');
    cfg.wave.headings=numeric_row(require_line(sections.PARA3,3,'PARA3'),'wave headings');
    if numel(cfg.wave.headings)~=cfg.wave.ndir
        error('CRESTU:HeadingCountMismatch','NDIR=%d but %d headings were supplied.', ...
            cfg.wave.ndir,numel(cfg.wave.headings));
    end

    v=numeric_row(require_line(sections.PARA4,1,'PARA4'),'PARA4'); require_count(v,3,'PARA4');
    cfg.water_depth=v(1); cfg.grav=v(2); cfg.rho=v(3);
    if cfg.water_depth<0 || cfg.grav<=0 || cfg.rho<=0
        error('CRESTU:InvalidEnvironment','Depth must be nonnegative; gravity and density must be positive.');
    end

    cfg.n_bodies=scalar_integer(require_line(sections.PARA5,1,'PARA5'),'NBODY');
    if cfg.n_bodies<1, error('CRESTU:InvalidBodyCount','NBODY must be at least one.'); end
    if numel(sections.PARA5)<cfg.n_bodies+1
        error('CRESTU:BodyConfigCount','PARA5 does not contain %d body records.',cfg.n_bodies);
    end
    cfg.bodies=repmat(struct('id',0,'mesh_file','','pos',zeros(1,3),'yaw',0),cfg.n_bodies,1);
    for b=1:cfg.n_bodies
        tokens=regexp(require_line(sections.PARA5,b+1,'PARA5'),'\s+','split');
        if numel(tokens)~=6, error('CRESTU:BodyRecordShape','Body record %d must contain 6 tokens.',b); end
        body_id=parse_scalar(tokens{1},sprintf('body %d ID',b)); mesh_name=tokens{2};
        if strcmpi(mesh_name,'auto')
            if cfg.n_bodies==1, mesh_name=sprintf('%s_body.bmf',cfg.case_name);
            else, mesh_name=sprintf('%s_body%d.bmf',cfg.case_name,body_id); end
        end
        cfg.bodies(b).id=body_id; cfg.bodies(b).mesh_file=resolve_path(cfg.config_dir,mesh_name);
        cfg.bodies(b).pos=cellfun(@(s)parse_scalar(s,'body position'),tokens(3:5));
        cfg.bodies(b).yaw=parse_scalar(tokens{6},'body yaw');
    end
    if numel(unique([cfg.bodies.id]))~=cfg.n_bodies
        error('CRESTU:DuplicateBodyId','Body IDs must be unique.');
    end

    if numel(sections.PARA6)<2*cfg.n_bodies
        error('CRESTU:MassConfigCount','PARA6 requires two lines per body.');
    end
    cfg.mass_props=repmat(struct('mass',0,'cg',zeros(1,3),'inertia',zeros(3)),cfg.n_bodies,1);
    for b=1:cfg.n_bodies
        m=numeric_row(sections.PARA6{2*b-1},sprintf('body %d mass/CG',b));
        q=numeric_row(sections.PARA6{2*b},sprintf('body %d inertia',b));
        require_count(m,4,'mass and CG'); require_count(q,6,'inertia tensor');
        cfg.mass_props(b).mass=m(1); cfg.mass_props(b).cg=m(2:4);
        cfg.mass_props(b).inertia=[q(1),-q(4),-q(5);-q(4),q(2),-q(6);-q(5),-q(6),q(3)];
        if cfg.mass_props(b).mass<=0 || any(eig(cfg.mass_props(b).inertia)<=0)
            error('CRESTU:InvalidMassProperties','Body %d mass/inertia must be positive definite.',b);
        end
    end

    cfg.calc_modes=numeric_row(require_line(sections.PARA7,1,'PARA7'),'PARA7 modes');
    if any(cfg.calc_modes<1|cfg.calc_modes>6|fix(cfg.calc_modes)~=cfg.calc_modes)
        error('CRESTU:InvalidModes','Radiation modes must be integers in [1,6].');
    end
    v=numeric_row(require_line(sections.PARA8,1,'PARA8'),'PARA8');
    require_count(v,2,'PARA8'); validate_flags(v,'PARA8'); cfg.isx=v(1); cfg.isy=v(2);
    % is x>=0 for ISX and y>=0 for ISY.  Polar-vector and axial-vector
    % parities follow reflection tensor rules.
    cfg.symmetry=struct();
    cfg.symmetry.multiplicity=2^(cfg.isx+cfg.isy);
    cfg.symmetry.area_scale=cfg.symmetry.multiplicity;
    cfg.symmetry.x_bounds=iff(cfg.isx==1,[0,Inf],[-Inf,Inf]);
    cfg.symmetry.y_bounds=iff(cfg.isy==1,[0,Inf],[-Inf,Inf]);
    cfg.symmetry.mode_parity_x=[-1,1,1,1,-1,-1]; % surge,sway,heave,roll,pitch,yaw
    cfg.symmetry.mode_parity_y=[1,-1,1,-1,1,-1];

    v=numeric_row(require_line(sections.PARA9,1,'PARA9'),'PARA9'); require_count(v,7,'PARA9');
    cfg.z_tol=v(1); cfg.fs=struct('nr_near',v(2),'nr_sponge',v(3),'r_inner',v(4), ...
        'r_outer',v(5),'sponge_ratio',v(6),'nz_farfield',v(7));
    cfg.fs.mu0=parse_scalar(require_line(sections.PARA10,1,'PARA10'),'PARA10');
    if cfg.z_tol<=0 || any([cfg.fs.nr_near,cfg.fs.nr_sponge,cfg.fs.nz_farfield]<1) || ...
            cfg.fs.r_inner<=0 || cfg.fs.r_outer<=cfg.fs.r_inner || ...
            cfg.fs.sponge_ratio<=0 || cfg.fs.mu0<0
        error('CRESTU:InvalidDomainParameters','PARA9/PARA10 parameters are invalid.');
    end

    cfg.files.fs=fullfile(cfg.config_dir,sprintf('%s_fs.bmf',cfg.case_name));
    cfg.files.seabed=fullfile(cfg.config_dir,sprintf('%s_seabed.bmf',cfg.case_name));
    cfg.files.farfield=fullfile(cfg.config_dir,sprintf('%s_farfield.bmf',cfg.case_name));
    cfg.files.potential_cache=fullfile(cfg.config_dir,sprintf('%s_PotCache.mat',cfg.case_name));
    cfg.files.results=fullfile(cfg.config_dir,sprintf('%s_Results.mat',cfg.case_name));
    fprintf('>>> Loaded %s\n',config_file);
    fprintf('    Case=%s, IPOTEN=%d, IFORCE=%d, IRAD=%d, IDIFF=%d, IDRIFT=%d\n',cfg.case_name, ...
        cfg.run.ipoten,cfg.run.iforce,cfg.run.irad,cfg.run.idiff,cfg.run.idrift);
    fprintf('    Bodies=%d, frequencies=%d, headings=%d, depth=%.3g m\n', ...
        cfg.n_bodies,cfg.freq.nfreq,cfg.wave.ndir,cfg.water_depth);
end

function sections=split_sections(lines)
% SPLIT_SECTIONS Execute the documented split_sections operation.
%
% Syntax:
%   sections=split_sections(lines)
%
% Inputs:
%   lines           : [documented value] Input required by the implemented function contract.
%
% Outputs:
%   sections        : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
    sections=struct(); current='';
    for k=1:numel(lines)
        token=regexp(lines{k},'^\s*(PARA\d+)\s*:','tokens','once','ignorecase');
        if ~isempty(token)
            current=upper(token{1});
            sections.(current)=cell(0,1);
        elseif isempty(current)
            error('CRESTU:ConfigPreamble','Data appears before a PARA section: %s',lines{k});
        else
            values=sections.(current);
            values{end+1,1}=lines{k}; %#ok<AGROW> section lengths are unknown while parsing
            sections.(current)=values;
        end
    end
end
function line=require_line(lines,idx,section)
% REQUIRE_LINE Execute the documented require_line operation.
%
% Syntax:
%   line=require_line(lines,idx,section)
%
% Inputs:
%   lines           : [documented value] Input required by the implemented function contract.
%   idx             : [integer scalar or array] Discrete count or index required by the algorithm.
%   section         : [documented value] Input required by the implemented function contract.
%
% Outputs:
%   line            : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
    if numel(lines)<idx, error('CRESTU:MissingConfigData','%s is missing line %d.',section,idx); end
    line=lines{idx};
end
function values=numeric_row(line,label)
% NUMERIC_ROW Execute the documented numeric_row operation.
%
% Syntax:
%   values=numeric_row(line,label)
%
% Inputs:
%   line            : [documented value] Input required by the implemented function contract.
%   label           : [documented value] Input required by the implemented function contract.
%
% Outputs:
%   values          : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
    tokens=regexp(strtrim(normalize_text(line)),'[,\s]+','split'); values=zeros(1,numel(tokens));
    for k=1:numel(tokens), values(k)=parse_scalar(tokens{k},label); end
end
function value=parse_scalar(token,label)
% PARSE_SCALAR Execute the documented parse_scalar operation.
%
% Syntax:
%   value=parse_scalar(token,label)
%
% Inputs:
%   token           : [documented value] Input required by the implemented function contract.
%   label           : [documented value] Input required by the implemented function contract.
%
% Outputs:
%   value           : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
    token=normalize_text(token); value=str2double(token);
    if ~isscalar(value)||~isfinite(value), error('CRESTU:InvalidNumber','Invalid token "%s" in %s.',token,label); end
end
function value=scalar_integer(line,label)
% SCALAR_INTEGER Execute the documented scalar_integer operation.
%
% Syntax:
%   value=scalar_integer(line,label)
%
% Inputs:
%   line            : [documented value] Input required by the implemented function contract.
%   label           : [documented value] Input required by the implemented function contract.
%
% Outputs:
%   value           : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
    values=numeric_row(line,label); require_count(values,1,label); value=values(1);
    if value<0||fix(value)~=value, error('CRESTU:InvalidInteger','%s must be a nonnegative integer.',label); end
end
function require_count(values,count,label)
% REQUIRE_COUNT Execute the documented require_count operation.
%
% Syntax:
%   require_count(values,count,label)
%
% Inputs:
%   values          : [documented value] Input required by the implemented function contract.
%   count           : [integer scalar or array] Discrete count or index required by the algorithm.
%   label           : [documented value] Input required by the implemented function contract.
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
    if numel(values)~=count
        error('CRESTU:ConfigArrayShape','%s requires %d values; got %d.',label,count,numel(values));
    end
end
function validate_flags(values,label)
% VALIDATE_FLAGS Execute the documented validate_flags operation.
%
% Syntax:
%   validate_flags(values,label)
%
% Inputs:
%   values          : [documented value] Input required by the implemented function contract.
%   label           : [documented value] Input required by the implemented function contract.
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
    if any(values~=0&values~=1), error('CRESTU:InvalidFlag','%s flags must be 0 or 1.',label); end
end
function out=normalize_text(in)
% NORMALIZE_TEXT Execute the documented normalize_text operation.
%
% Syntax:
%   out=normalize_text(in)
%
% Inputs:
%   in              : [documented value] Input required by the implemented function contract.
%
% Outputs:
%   out             : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
    out=strrep(in,char(8722),'-'); out=strrep(out,char(8208),'-');
    out=strrep(out,char(8209),'-'); out=strrep(out,char(8211),'-'); out=strrep(out,char(160),' ');
end
function out=absolute_path(in)
% ABSOLUTE_PATH Execute the documented absolute_path operation.
%
% Syntax:
%   out=absolute_path(in)
%
% Inputs:
%   in              : [documented value] Input required by the implemented function contract.
%
% Outputs:
%   out             : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
    if is_absolute(in), out=in; else, out=fullfile(pwd,in); end
end
function out=resolve_path(base_dir,in)
% RESOLVE_PATH Execute the documented resolve_path operation.
%
% Syntax:
%   out=resolve_path(base_dir,in)
%
% Inputs:
%   base_dir        : [documented value] Input required by the implemented function contract.
%   in              : [documented value] Input required by the implemented function contract.
%
% Outputs:
%   out             : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
    if is_absolute(in), out=in; else, out=fullfile(base_dir,in); end
end
function tf=is_absolute(in)
% IS_ABSOLUTE Execute the documented is_absolute operation.
%
% Syntax:
%   tf=is_absolute(in)
%
% Inputs:
%   in              : [documented value] Input required by the implemented function contract.
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
    tf=~isempty(regexp(in,'^[A-Za-z]:[\\/]','once'))||startsWith(in,filesep);
end
function out=iff(condition,a,b)
% IFF Execute the documented iff operation.
%
% Syntax:
%   out=iff(condition,a,b)
%
% Inputs:
%   condition       : [documented value] Input required by the implemented function contract.
%   a               : [documented value] Input required by the implemented function contract.
%   b               : [documented value] Input required by the implemented function contract.
%
% Outputs:
%   out             : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
    if condition, out=a; else, out=b; end
end
