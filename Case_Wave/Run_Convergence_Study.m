function study = Run_Convergence_Study(run_solver)
%RUN_CONVERGENCE_STUDY Three-grid hemisphere study against WAMIT IRR0/IRR3.
% Total-domain targets: approximately 500/1500/3500 panels.
% Full-domain meshes are used so symmetry acceleration cannot bias the grid-error assessment.
    if nargin<1, run_solver=true; end
    here=fileparts(mfilename('fullpath')); root=fileparts(here);
    addpath(fullfile(root,'1.Input'),fullfile(root,'2.Mesh'),fullfile(root,'3.HessSmith'), ...
        fullfile(root,'4.Potential'),fullfile(root,'5.Force'),fullfile(root,'6.MeanDriftLoads'),here);
    study_file=fullfile(here,'Convergence_Study.mat');
    if ~run_solver&&exist(study_file,'file')
        loaded=load(study_file,'study'); study=loaded.study; return;
    end
    levels=struct( ...
        'name',{'Coarse','Medium','Fine'},'N',{3,5,7}, ...
        'config',{'CRESTU_Convergence_Coarse.cfg','CRESTU_Convergence_Medium.cfg','CRESTU_Convergence_Fine.cfg'}, ...
        'body_file',{'hemi_D10_conv_coarse_body.bmf','hemi_D10_conv_medium_body.bmf','hemi_D10_conv_fine_body.bmf'});
    results=cell(3,1); runtimes=zeros(3,1);
    for level=1:3
        body_file=fullfile(here,levels(level).body_file); config_file=fullfile(here,levels(level).config);
        if run_solver
            generate_body_bmf(body_file,10,0,0,levels(level).N);
            timer=tic; results{level}=HydroMain(config_file); runtimes(level)=toc(timer);
        else
            cfg=read_config(config_file);
            if ~exist(cfg.files.results,'file')
                error('CRESTU:MissingConvergenceResult', ...
                    'Missing %s. Call Run_Convergence_Study(true).',cfg.files.results);
            end
            loaded=load(cfg.files.results,'results'); results{level}=loaded.results;
        end
    end
    wamit_root=fullfile(root,'..','WAMIT');
    reference=struct();
    reference.IRR0=load_reference(fullfile(wamit_root,'FullSphereIRR0'),results{1}.config);
    reference.IRR3=load_reference(fullfile(wamit_root,'FullSphereIRR3'),results{1}.config);
    plot_file=fullfile(here,'Convergence_Study.png');
    drift_plot_file=fullfile(here,'MeanDrift_Validation.png');
    make_first_order_plot(results,levels,reference,plot_file);
    make_drift_plot(results{3},reference,drift_plot_file);
    summary=make_summary(results,levels,runtimes);
    study=struct('levels',levels,'results',{results},'runtimes',runtimes, ...
        'reference',reference,'summary',summary,'plot_file',plot_file, ...
        'drift_plot_file',drift_plot_file);
    writetable(summary,fullfile(here,'Convergence_Study.csv'));
    save(study_file,'study','-v7.3');
end

function reference=load_reference(folder,cfg)
    reference.first_order=read_wamit_first_order(fullfile(folder,'hemisphere.1'),cfg.rho);
    reference.excitation=read_wamit_excitation(fullfile(folder,'hemisphere.2'),cfg.rho,cfg.grav,1,1);
    reference.rao=read_wamit_rao(fullfile(folder,'hemisphere.4'));
    reference.drift_momentum=read_wamit_mean_drift(fullfile(folder,'hemisphere.8'),cfg.rho,cfg.grav,1,1);
    reference.drift_pressure=read_wamit_mean_drift(fullfile(folder,'hemisphere.9'),cfg.rho,cfg.grav,1,1);
end

function make_first_order_plot(results,levels,reference,filename)
    figure_handle=figure('Visible','off','Color','w','Position',[100,100,1200,820]);
    layout=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    labels={levels.name}; colors=lines(3);
    axes_list=gobjects(4,1);
    for p=1:4, axes_list(p)=nexttile(layout); hold(axes_list(p),'on'); grid(axes_list(p),'on'); end
    for level=1:3
        r=results{level}; w=r.omegas;
        plot(axes_list(1),w,squeeze(r.added_mass(3,3,:)),'o-','Color',colors(level,:),'DisplayName',labels{level});
        plot(axes_list(2),w,squeeze(r.damping(3,3,:)),'o-','Color',colors(level,:),'DisplayName',labels{level});
        plot(axes_list(3),w,abs(squeeze(r.excitation(3,1,:))),'o-','Color',colors(level,:),'DisplayName',labels{level});
        plot(axes_list(4),w,squeeze(r.rao.amplitude(3,1,:)),'o-','Color',colors(level,:),'DisplayName',labels{level});
    end
    styles={'k--','k:'}; names={'WAMIT IRR0','WAMIT IRR3'}; refs={reference.IRR0,reference.IRR3};
    for q=1:2
        ref=refs{q}; [wr,order]=sort(ref.first_order.omegas);
        plot(axes_list(1),wr,squeeze(ref.first_order.added_mass(3,3,order)),styles{q},'DisplayName',names{q});
        plot(axes_list(2),wr,squeeze(ref.first_order.damping(3,3,order)),styles{q},'DisplayName',names{q});
        [we,order]=sort(ref.excitation.omegas);
        plot(axes_list(3),we,abs(squeeze(ref.excitation.force(3,1,order))), ...
            styles{q},'DisplayName',names{q});
        [wx,order]=sort(ref.rao.omegas);
        plot(axes_list(4),wx,squeeze(ref.rao.amplitude(3,1,order)), ...
            styles{q},'DisplayName',names{q});
    end
    ylabel(axes_list(1),'A_{33} [kg]'); ylabel(axes_list(2),'B_{33} [kg/s]');
    ylabel(axes_list(3),'|F_{z,exc}| [N/m]'); ylabel(axes_list(4),'|\xi_3|/A');
    for p=1:4, xlabel(axes_list(p),'\omega [rad/s]'); xlim(axes_list(p),[0.7,1.3]); end
    legend(axes_list(1),'Location','best'); title(layout,'CRESTU hemisphere grid convergence');
    exportgraphics(figure_handle,filename,'Resolution',180); close(figure_handle);
end

function make_drift_plot(result,reference,filename)
    if ~isfield(result,'drift')||~result.drift.enabled, return; end
    norm_scale=result.drift.normalization; w=result.omegas;
    near=squeeze(result.drift.near(1,1,:))/norm_scale;
    far=squeeze(result.drift.far(1,1,:))/norm_scale;
    fig=figure('Visible','off','Color','w'); hold on; grid on;
    plot(w,near,'o-','DisplayName','CRESTU near field');
    plot(w,far,'s-','DisplayName','CRESTU far field');
    refs={reference.IRR0,reference.IRR3}; styles={'k--','k:'}; names={'WAMIT IRR0 (.8/.9)','WAMIT IRR3 (.8/.9)'};
    for q=1:2
        direct=refs{q}.drift_pressure; momentum=refs{q}.drift_momentum;
        [wd,order]=sort(direct.omegas);
        plot(wd,squeeze(direct.Cd(1,1,order)),styles{q}, ...
            'DisplayName',[names{q},' pressure']);
        [wm,order]=sort(momentum.omegas);
        plot(wm,squeeze(momentum.Cd(1,1,order)),styles{q}, ...
            'LineWidth',1.5,'DisplayName',[names{q},' momentum']);
    end
    xlabel('\omega [rad/s]'); ylabel('C_{d,1}'); xlim([0.7,1.3]);
    title('Mean surge drift: direct pressure vs Kochin momentum'); legend('Location','best');
    exportgraphics(fig,filename,'Resolution',180); close(fig);
end

function summary=make_summary(results,levels,runtimes)
    rows=0; for q=1:3, rows=rows+numel(results{q}.omegas); end
    grid_name=strings(rows,1); omega=zeros(rows,1); computational=zeros(rows,1); equivalent=zeros(rows,1);
    A33=zeros(rows,1); B33=zeros(rows,1); Fz=zeros(rows,1); heave_rao=zeros(rows,1); runtime=zeros(rows,1); cursor=0;
    for q=1:3
        r=results{q}; n=numel(r.omegas); idx=cursor+(1:n); grid_name(idx)=levels(q).name; omega(idx)=r.omegas;
        computational(idx)=r.stats.total_dofs; equivalent(idx)=r.symmetry.full_equivalent_dofs;
        A33(idx)=squeeze(r.added_mass(3,3,:)); B33(idx)=squeeze(r.damping(3,3,:));
        Fz(idx)=abs(squeeze(r.excitation(3,1,:))); heave_rao(idx)=squeeze(r.rao.amplitude(3,1,:));
        runtime(idx)=runtimes(q); cursor=cursor+n;
    end
    summary=table(grid_name,omega,computational,equivalent,A33,B33,Fz,heave_rao,runtime);
end
