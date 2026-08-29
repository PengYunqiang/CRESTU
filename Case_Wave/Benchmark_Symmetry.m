function benchmark = Benchmark_Symmetry(run_solver)
% BENCHMARK_SYMMETRY Compare reduced-domain symmetry solutions with a full-domain solution.
%
% Syntax:
%   benchmark = Benchmark_Symmetry(run_solver)
%
% Description:
%   Compare reduced-domain symmetry solutions with a full-domain solution.
%   The implementation preserves the CRESTU solver and data-field contracts.
%
% Inputs:
%   run_solver - Flag that enables a new solver run, logical scalar [-].
%
% Outputs:
%   benchmark - Symmetry benchmark result structure [-].
%
% Governing Equations / Theory:
%   Uses the linear potential-flow, source-panel, or post-processing
%   formulation stated in the CRESTU theory and technical manual.
%
% References:
%   - CRESTU theory and technical manual.
%   - Standard boundary-element and marine hydrodynamics references.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)
%BENCHMARK_SYMMETRY Exact-grid full, half and quarter symmetry comparison.
% A quarter mesh is reflected to construct the half/full reference meshes;
% therefore all three cases have identical panels after image expansion.
% The full-domain baseline is reflected from the same quarter mesh so that
% independently generated panel layouts cannot contaminate the comparison.
%% Stage 1: Initialize inputs and dependencies

    if nargin < 1
        run_solver = true;
    end

%% Stage 2: Run the core calculation
    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);
    addpath(fullfile(root,'1.Input'),fullfile(root,'2.Mesh'),fullfile(root,'3.HessSmith'), ...
        fullfile(root,'4.Potential'),fullfile(root,'5.Force'),fullfile(root,'6.MeanDriftLoads'),here);
    output_file = fullfile(here,'Symmetry_Benchmark.mat');
    if ~run_solver
        if ~exist(output_file,'file')
            error('CRESTU:MissingBenchmarkResult','Call Benchmark_Symmetry(true) first.');
        end
        loaded = load(output_file,'benchmark');
        benchmark = loaded.benchmark;
        return;
    end
    quarter = build_bmf_domain(fullfile(here,'CRESTU_SymQuarter.cfg'));
    domains = {expand_domain(quarter,1,1,0,0),expand_domain(quarter,0,1,1,0),quarter};
    labels = {'Full';'Half-X';'Quarter-XY'};
    cases = cell(3,1);
    runtime = zeros(3,1);
    for k = 1:3
        timer = tic;
        cases{k} = solve_domain_once(domains{k});
        runtime(k) = toc(timer);
    end
    unknowns = cellfun(@(r)r.stats.total_dofs,cases);
    equivalent = cellfun(@(r)r.symmetry.full_equivalent_dofs,cases);
    A33 = cellfun(@(r)r.added_mass(3,3),cases);
    B33 = cellfun(@(r)r.damping(3,3),cases);
    F3 = cellfun(@(r)r.excitation(3),cases);
    speedup = runtime(1) ./ max(runtime,eps);
    relA = abs(A33 - A33(1)) / max(abs(A33(1)),eps);
    relB = abs(B33 - B33(1)) / max(abs(B33(1)),eps);
    relF = abs(F3 - F3(1)) / max(abs(F3(1)),eps);
    summary = table(labels,unknowns,equivalent,runtime,speedup,A33,B33,F3,relA,relB,relF, ...
'VariableNames',{'domain','unknowns','full_equivalent_panels','runtime_s','speedup', ...
'A33','B33','F3','A33_rel_difference','B33_rel_difference','F3_rel_difference'});
    benchmark = struct('summary',summary,'cases',{cases}, ...
'note','Exact reflected-grid benchmark; runtime includes assembly and dense solves.');
    writetable(summary,fullfile(here,'Symmetry_Benchmark.csv'));
    save(output_file,'benchmark');
    disp(summary);
end

function expanded = expand_domain(base,reflect_x,reflect_y,target_isx,target_isy)
% EXPAND_DOMAIN Expand a reduced panel domain by the requested reflections.
%
% Syntax:
%   expanded=expand_domain(base,reflect_x,reflect_y,target_isx,target_isy)
%
% Description:
%   Expand a reduced panel domain by the requested reflections.
%   The implementation preserves the CRESTU solver and data-field contracts.
%
% Inputs:
%   base - Base boundary-domain structure [-].
%   reflect_x - Flag that expands the x reflection, logical scalar [-].
%   reflect_y - Flag that expands the y reflection, logical scalar [-].
%   target_isx - Target x-symmetry flag, logical scalar [-].
%   target_isy - Target y-symmetry flag, logical scalar [-].
%
% Outputs:
%   expanded - Expanded boundary-domain structure [-].
%
% Governing Equations / Theory:
%   Uses the linear potential-flow, source-panel, or post-processing
%   formulation stated in the CRESTU theory and technical manual.
%
% References:
%   - CRESTU theory and technical manual.
%   - Standard boundary-element and marine hydrodynamics references.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)
%% Stage 1: Initialize inputs and dependencies

    expanded = base;

%% Stage 2: Run the core calculation

    for bodyIndex = 1:numel(base.body_list)
        expanded.body_list{bodyIndex} = expand_mesh_by_symmetry(base.body_list{bodyIndex},reflect_x,reflect_y);
    end
    names = {'fs','seabed','farfield'};
    for k = 1:numel(names)
        name = names{k};
        if ~isempty(base.(name))
            expanded.(name) = expand_mesh_by_symmetry(base.(name),reflect_x,reflect_y);
        end
    end
    expanded.cfg.isx = target_isx;
    expanded.cfg.isy = target_isy;
    expanded.cfg.symmetry.multiplicity = 2^(target_isx + target_isy);
    expanded.cfg.symmetry.area_scale = expanded.cfg.symmetry.multiplicity;
    expanded.stats = domain_stats(expanded);
    expanded.geometry = merge_domain_geometry(expanded);
end

function stats = domain_stats(domain)
% DOMAIN_STATS Collect panel counts and system dimensions for one domain.
%
% Syntax:
%   stats=domain_stats(domain)
%
% Description:
%   Collect panel counts and system dimensions for one domain.
%   The implementation preserves the CRESTU solver and data-field contracts.
%
% Inputs:
%   domain - Prepared boundary-domain structure [-].
%
% Outputs:
%   stats - Panel-count and system-size summary structure [-].
%
% Governing Equations / Theory:
%   Uses the linear potential-flow, source-panel, or post-processing
%   formulation stated in the CRESTU theory and technical manual.
%
% References:
%   - CRESTU theory and technical manual.
%   - Standard boundary-element and marine hydrodynamics references.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)
%% Stage 1: Initialize inputs and dependencies

    stats = struct('total_body_panels',sum(cellfun(@(m)m.n_panels,domain.body_list)), ...
'fs_panels',panel_count(domain.fs),'seabed_panels',panel_count(domain.seabed), ...
'farfield_panels',panel_count(domain.farfield));
    stats.total_dofs = stats.total_body_panels + stats.fs_panels + stats.seabed_panels + stats.farfield_panels;
end

function panelCount = panel_count(mesh)
% PANEL_COUNT Return the number of panels in an optional mesh.
%
% Syntax:
%   panelCount=panel_count(mesh)
%
% Description:
%   Return the number of panels in an optional mesh.
%   The implementation preserves the CRESTU solver and data-field contracts.
%
% Inputs:
%   mesh - Panel-mesh structure or empty value [-].
%
% Outputs:
%   panelCount - Panel count, nonnegative integer [-].
%
% Governing Equations / Theory:
%   Uses the linear potential-flow, source-panel, or post-processing
%   formulation stated in the CRESTU theory and technical manual.
%
% References:
%   - CRESTU theory and technical manual.
%   - Standard boundary-element and marine hydrodynamics references.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)
%% Stage 1: Initialize inputs and dependencies

    if isempty(mesh)
        panelCount = 0;
    else
        panelCount = mesh.n_panels;
    end
end

function result = solve_domain_once(domain)
% SOLVE_DOMAIN_ONCE Solve one prepared radiation and diffraction domain.
%
% Syntax:
%   result=solve_domain_once(domain)
%
% Description:
%   Solve one prepared radiation and diffraction domain.
%   The implementation preserves the CRESTU solver and data-field contracts.
%
% Inputs:
%   domain - Prepared boundary-domain structure [-].
%
% Outputs:
%   result - Computed hydrodynamic result structure [-].
%
% Governing Equations / Theory:
%   Uses the linear potential-flow, source-panel, or post-processing
%   formulation stated in the CRESTU theory and technical manual.
%
% References:
%   - CRESTU theory and technical manual.
%   - Standard boundary-element and marine hydrodynamics references.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)
%% Stage 1: Initialize inputs and dependencies

    cfg = domain.cfg;

%% Stage 2: Run the core calculation

    geom = domain.geometry;
    omega = cfg.freq.omegas(1);
    heading = cfg.wave.headings(1);
    nb = domain.stats.total_body_panels;
    ndof = 6 * cfg.n_bodies;
    centers = geom.centers(1:nb,:);
    normals = geom.normals(1:nb,:);
    areas = geom.areas(1:nb);
    nj = compute_generalized_normals(centers,normals,domain.body_list);
    mode_parity = get_mode_parities(cfg.n_bodies,cfg.isx,cfg.isy);
    group_parity = unique(mode_parity,'rows','stable');
    ng = size(group_parity,1);
    phi_rad = complex(zeros(nb,ndof));
    force = complex(zeros(ndof,1));
    assembly = cell(ng,1);
    for groupIndex = 1:ng
        parity = group_parity(groupIndex,:);
        modes = find(all(mode_parity == parity,2));
        [K,S,assembly{groupIndex}] = assemble_rankine_matrix(geom.total_panels,geom.centers,geom.normals, ...
            geom.vertices,domain.stats,omega,cfg.grav,cfg,parity);
        [phi,factors] = solve_complex_system(K,S * (1i * omega * nj(:,modes)),[]);
        phi_rad(:,modes) = phi(1:nb,:);
        [phi_I,dphi_I] = decompose_incident_wave_symmetry(centers,normals,omega,cfg.grav, ...
            cfg.water_depth,heading,1,cfg.isx,cfg.isy,parity);
        phi_D = solve_complex_system(K,S * dphi_I,factors);
        component_force = compute_wave_excitation(phi_I,phi_D(1:nb),nj,areas,omega,cfg.rho);
        weights = symmetry_force_weights(mode_parity,parity,cfg.isx,cfg.isy);
        force = force + component_force .* weights;
    end
    symmetry = struct('isx',cfg.isx,'isy',cfg.isy,'multiplicity',2^(cfg.isx + cfg.isy), ...
'mode_parity',mode_parity,'full_equivalent_dofs',geom.total_panels * 2^(cfg.isx + cfg.isy));
    [A,pressure_damping,diagnostics] = compute_hydrodynamic_coeffs(phi_rad,nj,areas,omega,cfg.rho,symmetry);
    [B,energy_diagnostics] = compute_radiation_damping_energy(phi_rad,nj,centers,normals, ...
        areas,omega,cfg,mode_parity,72);
    diagnostics.pressure_damping = pressure_damping;
    diagnostics.energy_diagnostics = energy_diagnostics;
    result = struct('stats',domain.stats,'symmetry',symmetry,'added_mass',A,'damping',B, ...
'excitation',force,'diagnostics',diagnostics,'assembly_info',{assembly});
end
