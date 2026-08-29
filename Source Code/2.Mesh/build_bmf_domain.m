function domain = build_bmf_domain(config_file)
% BUILD_BMF_DOMAIN Build bmf domain for the CRESTU hydrodynamic workflow.
%
% Syntax:
%   domain = build_bmf_domain(config_file)
%
% Description:
%   Prepares boundary-panel geometry for the Rankine solver.
%   Global coordinates and panel-normal signs are preserved.
%
% Inputs:
%   config_file        - [character vector or string scalar] CRESTU configuration-file path.
%
% Outputs:
%   domain             - [struct] Complete boundary-element domain and component metadata in SI units.
%
% Governing Equations / Theory:
%   Planar panel geometry, polygon moments, reflection transformations, and mesh-topology relations as applicable.
%
% References:
%   - Hess and Smith (1964); CRESTU BMF mesh-format and geometry conventions.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    if nargin < 1 || isempty(config_file)
        config_file ='CRESTU.cfg';
    end
    cfg = read_config(config_file);
    [cfg, cfg.fs.tuning_report] = tune_sponge_layer(cfg, cfg.freq.omegas);
    body_list = cell(cfg.n_bodies, 1);
    waterlines = cell(cfg.n_bodies, 1);
    total_body_panels = 0;
    for bodyIndex = 1:cfg.n_bodies
        bc = cfg.bodies(bodyIndex);
        if cfg.isx && abs(bc.pos(1)) > cfg.z_tol
            error('CRESTU:SymmetryBodyMapping','ISX=1 currently requires each body reference position on x=0.');
        end
        if cfg.isy && abs(bc.pos(2)) > cfg.z_tol
            error('CRESTU:SymmetryBodyMapping','ISY=1 currently requires each body reference position on y=0.');
        end
        if ~exist(bc.mesh_file,'file')
            generate_body_bmf(bc.mesh_file, 10, cfg.isx, cfg.isy, 5);
        end
        mesh = transform_body_mesh(read_bmf(bc.mesh_file), bc);
        mesh = reduce_mesh_by_symmetry(mesh, cfg.isx, cfg.isy, cfg.z_tol);
        mesh.cg = cfg.mass_props(bodyIndex).cg;
        mesh.body_id = bc.id;
        body_list{bodyIndex} = mesh;
        waterlines{bodyIndex} = extract_waterline(mesh, cfg.z_tol);
        total_body_panels = total_body_panels + mesh.n_panels;
    end
    if cfg.n_bodies == 1
        waterline = waterlines{1};
        mesh_fs = generate_free_surface_bmf(cfg.files.fs, waterline, cfg);
    else
        waterline = waterlines;
        mesh_fs = generate_multibody_free_surface_bmf(cfg.files.fs, waterlines, cfg);
    end
    has_seabed = cfg.water_depth > 0;
    mesh_seabed = [];
    mesh_farfield = [];
    reference_waterline = waterlines{1};
    if has_seabed
        if (cfg.isx || cfg.isy) && cfg.n_bodies == 1
            full_wl = complete_waterline_by_symmetry(reference_waterline, cfg.isx, cfg.isy);
            cfg_full = cfg;
            cfg_full.isx = 0;
            cfg_full.isy = 0;
            mesh_seabed = generate_reduced_seabed_mesh(cfg, full_wl.n_pts);
            mesh_farfield = generate_farfield_mesh(full_wl, cfg_full, cfg.fs.nz_farfield);
            mesh_farfield = reduce_mesh_by_symmetry(mesh_farfield, cfg.isx, cfg.isy, cfg.z_tol);
        elseif cfg.n_bodies == 1
            mesh_seabed = generate_seabed_mesh(reference_waterline, cfg, mesh_fs);
            mesh_farfield = generate_farfield_mesh(reference_waterline, cfg, cfg.fs.nz_farfield);
        else
            mesh_seabed = generate_seabed_disk_mesh(cfg, mesh_fs);
            mesh_farfield = generate_farfield_mesh(reference_waterline, cfg, cfg.fs.nz_farfield);
        end
        write_bmf(cfg.files.seabed, mesh_seabed);
        write_bmf(cfg.files.farfield, mesh_farfield);
    end
    stats = struct('total_body_panels', total_body_panels,'fs_panels', mesh_fs.n_panels, ...
'seabed_panels', panel_count(mesh_seabed),'farfield_panels', panel_count(mesh_farfield));
    stats.total_dofs = stats.total_body_panels + stats.fs_panels + stats.seabed_panels + stats.farfield_panels;
    domain = struct('cfg', cfg,'body_list', {body_list},'fs', mesh_fs,'seabed', mesh_seabed, ...
'farfield', mesh_farfield,'waterline', {waterline},'stats', stats);
    domain.geometry = merge_domain_geometry(domain);
    fprintf('[INFO] Domain: %d bodies / %d body panels / %d total panels.\n', ...
        cfg.n_bodies, total_body_panels, stats.total_dofs);
end

function panelCount = panel_count(mesh)
% PANEL_COUNT Return the number of panels stored in a mesh structure.
%
% Syntax:
%   panelCount = panel_count(mesh)
%
% Description:
%   Prepares boundary-panel geometry for the Rankine solver.
%   Global coordinates and panel-normal signs are preserved.
%
% Inputs:
%   mesh               - [struct] Boundary-panel mesh with Cartesian geometry in SI units.
%
% Outputs:
%   panelCount                  - [scalar] Number of panels, degrees of freedom, or rows requested by the calling function, dimensionless.
%
% Governing Equations / Theory:
%   Planar panel geometry, polygon moments, reflection transformations, and mesh-topology relations as applicable.
%
% References:
%   - Hess and Smith (1964); CRESTU BMF mesh-format and geometry conventions.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    if isempty(mesh)
        panelCount = 0;
    else
        panelCount = mesh.n_panels;
    end
end
