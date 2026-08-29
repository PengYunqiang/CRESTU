function report = Test_MultiBody_Wave(run_full_domain)
% TEST_MULTIBODY_WAVE Compare reduced and full-domain multibody wave solutions.
%
% Syntax:
%   report = Test_MultiBody_Wave(run_full_domain)
%
% Description:
%   Compare reduced and full-domain multibody wave solutions.
%   The implementation preserves the CRESTU solver and data-field contracts.
%
% Inputs:
%   run_full_domain - Flag that enables the full-domain comparison, logical scalar [-].
%
% Outputs:
%   report - Multibody comparison report structure [-].
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
%TEST_MULTIBODY_WAVE 12-DOF two-sphere BEM/force/RAO regression.
% Default uses two coarse body meshes in an unbounded-fluid BEM smoke test.
% Pass true to run the configured free-surface domain through HydroMain.
%% Stage 1: Initialize inputs and dependencies

    if nargin < 1
        run_full_domain = false;
    end

%% Stage 2: Run the core calculation
    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);
    addpath(fullfile(root,'1.Input'),fullfile(root,'2.Mesh'),fullfile(root,'3.HessSmith'), ...
        fullfile(root,'4.Potential'),fullfile(root,'5.Force'),here);
    config_file = fullfile(here,'CRESTU_MultiBody.cfg');
    if run_full_domain
        results = HydroMain(config_file);
        assert(isequal(size(results.added_mass,1),12),'Expected a 12-DOF result.');
        report = struct('mode','full_domain','passed',true,'results',results);
        return;
    end
    cfg = read_config(config_file);
    coarse_file = fullfile(here,'two_sphere_smoke_body.bmf');
    generate_body_bmf(coarse_file,10,0,0,2);
    bodies = cell(2,1);
    for bodyIndex = 1:2
        local = read_bmf(coarse_file);
        bodies{bodyIndex} = transform_body_mesh(local,cfg.bodies(bodyIndex));
        bodies{bodyIndex}.cg = cfg.mass_props(bodyIndex).cg;
        bodies{bodyIndex}.body_id = bodyIndex;
    end
    domain = struct('body_list',{bodies},'fs',[],'seabed',[],'farfield',[]);
    geom = merge_domain_geometry(domain);
    stats = struct('total_body_panels',geom.body_panels,'fs_panels',0,'seabed_panels',0, ...
'farfield_panels',0,'total_dofs',geom.total_panels);
    nj = compute_generalized_normals(geom.centers,geom.normals,bodies);
    assert(isequal(size(nj),[geom.body_panels,12]),'Generalized-normal matrix is not Npanel-by-12.');
    omega = cfg.freq.omegas(1);
    [system,S] = assemble_rankine_matrix(geom.total_panels,geom.centers,geom.normals,geom.vertices, ...
        stats,omega,cfg.grav,cfg);
    phi = solve_complex_system(system,S * (1i * omega * nj));
    [A,B,diag_info] = compute_hydrodynamic_coeffs(phi,nj,geom.areas,omega,cfg.rho);
    phi_I = compute_incident_wave(geom.centers,geom.normals,omega,cfg.grav,0,0,1);
    [F,p] = compute_wave_excitation(phi_I,zeros(size(phi_I)),nj,geom.areas,omega,cfg.rho);
    [C,hydro] = compute_hydrostatic_matrix(bodies,cfg);
    rao = solve_rao(omega,A,B,C,reshape(F,12,1,1),cfg);
    assert(isequal(size(A),[12,12]) && isequal(size(B),[12,12]) && isequal(size(C),[12,12]));
    assert(all(isfinite([A(:);B(:);C(:);F(:);rao.complex(:);p(:)])),'Nonfinite multi-body result.');
    report = struct('mode','body_only_dimension_smoke','passed',true,'n_panels',geom.total_panels, ...
'generalized_normal_size',size(nj),'matrix_size',size(system),'added_mass',A, ...
'damping',B,'excitation',F,'hydrostatic',C,'hydrostatic_details',hydro, ...
'rao',rao,'diagnostics',diag_info);
    fprintf('[OK] Multi-body regression passed: %d panels, BEM %dx%d, 12 radiation modes.\n', ...
        geom.total_panels,size(system,1),size(system,2));
end
