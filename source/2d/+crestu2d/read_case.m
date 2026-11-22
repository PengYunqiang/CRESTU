function cfg = read_case(input)
%READ_CASE Explicit SI JSON input, or the same resolved structure in memory.
    if ischar(input) || isstring(input)
        assert(text_scalar(input),'CRESTU2D:Input','Input filename must be nonempty scalar text.');
        filename=char(input);
        assert(isfile(filename),'CRESTU2D:Input','Input file missing: %s',filename);
        [~,info]=fileattrib(filename); filename=info.Name;
        cfg=jsondecode(fileread(filename)); base=fileparts(filename);
        cfg.input_file=filename;
    else
        assert(isstruct(input) && isscalar(input),'CRESTU2D:Input','Expected JSON filename or scalar structure.');
        cfg=input; base=pwd;
    end
    assert(strcmp(cfg.schema,'CRESTU-2D-1'),'CRESTU2D:Schema','Unsupported input schema.');
    assert(strcmp(cfg.problem,'radiation_diffraction') || strcmp(cfg.problem,'radiation'), ...
        'CRESTU2D:Problem','Unsupported problem.');
    assert(strcmp(cfg.coordinates,'transverse_yz_z_up'),'CRESTU2D:Coordinates','Only transverse (y,z), z up is supported.');
    assert(strcmp(cfg.phase,'exp(+i*omega*t)'),'CRESTU2D:Phase','Expected exp(+i*omega*t).');
    assert(real_scalar(cfg.fluid.rho_kg_m3) && real_scalar(cfg.fluid.g_m_s2) && ...
        min([cfg.fluid.rho_kg_m3 cfg.fluid.g_m_s2])>0,'CRESTU2D:Fluid','rho and g must be positive.');
    assert(any(strcmp(cfg.fluid.domain,{'deep_water','unbounded','finite_depth'})), ...
        'CRESTU2D:Domain','Supported domains: deep_water, unbounded and finite_depth.');
    finiteDepth=strcmp(cfg.fluid.domain,'finite_depth');
    if finiteDepth
        assert(real_scalar(cfg.fluid.depth_m) && cfg.fluid.depth_m>0, ...
            'CRESTU2D:Depth','finite_depth requires positive numeric fluid.depth_m in metres.');
    else
        assert(text_scalar(cfg.fluid.depth_m) && strcmp(cfg.fluid.depth_m,'infinite'), ...
            'CRESTU2D:Depth','Deep/unbounded models require depth_m="infinite".');
    end
    if strcmp(cfg.fluid.domain,'unbounded')
        assert(strcmp(cfg.problem,'radiation'),'CRESTU2D:Domain','Unbounded fluid supports radiation only.');
    end
    assert(isnumeric(cfg.omega_rad_s) && isreal(cfg.omega_rad_s) && isvector(cfg.omega_rad_s) && ...
        ~isempty(cfg.omega_rad_s) && all(isfinite(cfg.omega_rad_s)) && all(cfg.omega_rad_s>0), ...
        'CRESTU2D:Frequency','Supply strictly positive angular frequencies in rad/s.');
    cfg.omega_rad_s=cfg.omega_rad_s(:).';
    assert(isnumeric(cfg.reference_yz_m) && isreal(cfg.reference_yz_m) && ...
        numel(cfg.reference_yz_m)==2 && all(isfinite(cfg.reference_yz_m)), ...
        'CRESTU2D:Reference','reference_yz_m must contain finite [y,z].');
    cfg.reference_yz_m=cfg.reference_yz_m(:).';
    assert(real_scalar(cfg.waterline_z_m),'CRESTU2D:Waterline','Invalid waterline.');
    assert(real_scalar(cfg.normalization_length_m) && cfg.normalization_length_m>0, ...
        'CRESTU2D:Scale','Positive normalization length required.');
    dofs=cellstr(string(cfg.dofs)); dofs=dofs(:).';
    [ok,ids]=ismember(dofs,{'sway','heave','roll'});
    assert(all(ok) && numel(unique(ids))==numel(ids) && ~isempty(ids),'CRESTU2D:DOF','Invalid or repeated DOF.');
    cfg.dofs=dofs; cfg.dof_indices=ids(:).';
    if finiteDepth
        assert(strcmp(cfg.solver.method,'direct_laplace_modal_dtn'), ...
            'CRESTU2D:Method','Finite depth requires direct_laplace_modal_dtn; the deep Green kernel is not a bottom model.');
        for name={'side_extent_m','surface_step_m','bottom_probe_step_m'}
            q=cfg.solver.(name{1});assert(real_scalar(q) && q>0,'CRESTU2D:Discretization','Finite-domain distances must be positive.');
        end
        q=cfg.solver.side_panels;
        assert(real_scalar(q) && q==fix(q) && q>=2,'CRESTU2D:Discretization','side_panels must be an integer >=2.');
    else
        for name={'wave_order','force_order'}
            q=cfg.solver.(name{1});
            assert(real_scalar(q) && q==fix(q) && q>=2,'CRESTU2D:Quadrature','Orders must be integers >=2.');
        end
        assert(strcmp(cfg.solver.method,'constant_source_exact_log'), ...
            'CRESTU2D:Method','Unsupported solver method.');
    end
    if ~isfield(cfg.solver,'augmentation'),cfg.solver.augmentation='none';end
    assert(~finiteDepth || strcmp(cfg.solver.augmentation,'none'), ...
        'CRESTU2D:Augmentation','Finite direct-domain BEM has its own representation; deep-water augmentation must be none.');
    assert(text_scalar(cfg.solver.augmentation) && any(strcmp(cfg.solver.augmentation,{'none','one_source_constraint','rigid_lid'})), ...
        'CRESTU2D:Augmentation','augmentation must be none, one_source_constraint or rigid_lid.');
    if strcmp(cfg.solver.augmentation,'one_source_constraint')
        assert(strcmp(cfg.fluid.domain,'deep_water'),'CRESTU2D:Augmentation','Waterline-source augmentation requires deep water.');
        assert(real_scalar(cfg.solver.auxiliary_source_y_m),'CRESTU2D:Augmentation','Supply one finite waterline source y coordinate.');
        p=cfg.solver.interior_constraint_yz_m;
        assert(isnumeric(p) && isreal(p) && numel(p)==2 && all(isfinite(p)), ...
            'CRESTU2D:Augmentation','Supply one finite [y,z] interior representation constraint.');
        cfg.solver.interior_constraint_yz_m=p(:).';
    end
    if strcmp(cfg.solver.augmentation,'rigid_lid')
        assert(strcmp(cfg.fluid.domain,'deep_water'),'CRESTU2D:Augmentation','Internal lid augmentation requires deep water.');
        assert(isfield(cfg.solver,'lid_panels') && real_scalar(cfg.solver.lid_panels) && ...
            cfg.solver.lid_panels==fix(cfg.solver.lid_panels) && cfg.solver.lid_panels>=2, ...
            'CRESTU2D:Augmentation','Supply an explicit integer lid_panels >=2.');
        if ~isfield(cfg.solver,'lid_quadrature_order'),cfg.solver.lid_quadrature_order=max(32,cfg.solver.wave_order);end
        assert(real_scalar(cfg.solver.lid_quadrature_order) && cfg.solver.lid_quadrature_order==fix(cfg.solver.lid_quadrature_order) && ...
            cfg.solver.lid_quadrature_order>=2,'CRESTU2D:Quadrature','Invalid lid_quadrature_order.');
    end
    if strcmp(cfg.problem,'radiation_diffraction')
        assert(real_scalar(cfg.incident.direction_y) && ismember(cfg.incident.direction_y,[-1 1]), ...
            'CRESTU2D:Wave','direction_y must be one real scalar, -1 or +1.');
        assert(real_scalar(cfg.incident.amplitude_m) && cfg.incident.amplitude_m>0, ...
            'CRESTU2D:Wave','Positive incident surface-elevation amplitude required.');
    end
    assert(text_scalar(cfg.name),'CRESTU2D:Input','name must be nonempty scalar text.');
    cfg.name=char(cfg.name);
    assert(flag_scalar(cfg.geometry.check_self_intersection),'CRESTU2D:Geometry', ...
        'check_self_intersection must be scalar true/false.');
    assert(flag_scalar(cfg.output.save) && text_scalar(cfg.output.directory),'CRESTU2D:Output', ...
        'output.save must be scalar true/false and output.directory nonempty scalar text.');
    cfg.output.save=logical(cfg.output.save);
    cfg.output.directory=char(cfg.output.directory);
    if ~isfield(cfg.geometry,'nodes_yz_m')
        assert(text_scalar(cfg.geometry.nodes_file),'CRESTU2D:Geometry','nodes_file must be nonempty scalar text.');
        path=resolve(base,cfg.geometry.nodes_file);
        assert(isfile(path),'CRESTU2D:Geometry','Section nodes file missing: %s',path);
        cfg.geometry.nodes_yz_m=readmatrix(path);
        cfg.geometry.nodes_file=path;
    end
    cfg.output.directory=resolve(base,cfg.output.directory);
end
function ok=real_scalar(value)
    ok=isnumeric(value) && isreal(value) && isscalar(value) && isfinite(value);
end
function ok=flag_scalar(value)
    ok=(islogical(value) || isnumeric(value)) && isreal(value) && isscalar(value) && any(value==[0 1]);
end
function ok=text_scalar(value)
    ok=(ischar(value) && isrow(value) && ~isempty(value)) || ...
        (isstring(value) && isscalar(value) && ~ismissing(value) && strlength(value)>0);
end
function path=resolve(base,path)
    path=char(path);
    if ~(startsWith(path,'/') || startsWith(path,'\\') || (numel(path)>=2 && path(2)==':'))
        path=fullfile(base,path);
    end
end
