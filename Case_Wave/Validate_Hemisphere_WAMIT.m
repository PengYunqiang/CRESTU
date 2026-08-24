function validation = Validate_Hemisphere_WAMIT(run_solver,config_file)
%VALIDATE_HEMISPHERE_WAMIT Compare saved/new hemi_D10 results with WAMIT .1.
    if nargin<1, run_solver=false; end
    here=fileparts(mfilename('fullpath')); root=fileparts(here);
    addpath(fullfile(root,'1.Input'),fullfile(root,'2.Mesh'), ...
        fullfile(root,'4.Potential'),fullfile(root,'5.Force'),here);
    if nargin<2||isempty(config_file)
        if run_solver, config_file=fullfile(here,'CRESTU_Validation.cfg');
        else, config_file=fullfile(here,'CRESTU.cfg'); end
    end
    cfg=read_config(config_file);
    if run_solver
        results=HydroMain(cfg.config_file);
    elseif exist(cfg.files.results,'file')
        data=load(cfg.files.results,'results'); results=data.results;
    else
        results=[];
    end
    wamit_file=fullfile(root,'..','WAMIT','FullSphereIRR0','hemisphere.1');
    reference=read_wamit_first_order(wamit_file,cfg.rho);
    validation=struct('reference',reference,'solver_results_available',~isempty(results), ...
        'comparison',[],'sanity',struct());
    if isempty(results)
        warning('CRESTU:NoSolverResults','No hemi_D10 result cache exists; WAMIT reference ingestion only.'); return;
    end
    validation.comparison=compare_wamit_results(results,reference);
    A33=squeeze(results.added_mass(3,3,:)); B33=squeeze(results.damping(3,3,:));
    validation.sanity=struct('A33_positive',all(A33>0),'B33_positive',all(B33>=0), ...
        'A_symmetry_max',max(arrayfun(@(k)norm(results.added_mass(:,:,k)-results.added_mass(:,:,k).','fro')/ ...
            max(norm(results.added_mass(:,:,k),'fro'),eps),1:numel(results.omegas))), ...
        'B_symmetry_max',max(arrayfun(@(k)norm(results.damping(:,:,k)-results.damping(:,:,k).','fro')/ ...
            max(norm(results.damping(:,:,k),'fro'),eps),1:numel(results.omegas))), ...
        'low_frequency_B33',B33(1),'high_frequency_B33',B33(end));
end
