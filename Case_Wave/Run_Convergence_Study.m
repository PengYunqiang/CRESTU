function study = Run_Convergence_Study(run_solver)
% RUN_CONVERGENCE_STUDY Run and report the three-grid hemisphere convergence study.
%
% Syntax:
%   study = Run_Convergence_Study(run_solver)
%
% Description:
%   Run and report the three-grid hemisphere convergence study.
%   The implementation preserves the CRESTU solver and data-field contracts.
%
% Inputs:
%   run_solver - Flag that enables a new solver run, logical scalar [-].
%
% Outputs:
%   study - Convergence-study result structure [-].
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
%RUN_CONVERGENCE_STUDY Three-grid hemisphere study against WAMIT IRR0/IRR3.
% Total-domain targets: approximately 500/1500/3500 panels.
% Full-domain meshes are used so symmetry acceleration cannot bias the grid-error assessment.
%% Stage 1: Initialize inputs and dependencies

    if nargin < 1
        run_solver = true;
    end

%% Stage 2: Run the core calculation
    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);
    addpath(fullfile(root,'1.Input'),fullfile(root,'2.Mesh'),fullfile(root,'3.HessSmith'), ...
        fullfile(root,'4.Potential'),fullfile(root,'5.Force'),fullfile(root,'6.MeanDriftLoads'),here);
    study_file = fullfile(here,'Convergence_Study.mat');
    if ~run_solver && exist(study_file,'file')
        loaded = load(study_file,'study');
        study = loaded.study;
        return;
    end
    levels = struct( ...
'name',{'Coarse','Medium','Fine'},'N',{3,5,7}, ...
'config',{'CRESTU_Convergence_Coarse.cfg','CRESTU_Convergence_Medium.cfg','CRESTU_Convergence_Fine.cfg'}, ...
'body_file',{'hemi_D10_conv_coarse_body.bmf','hemi_D10_conv_medium_body.bmf','hemi_D10_conv_fine_body.bmf'});
    results = cell(3,1);
    runtimes = zeros(3,1);
    for level = 1:3
        body_file = fullfile(here,levels(level).body_file);
        config_file = fullfile(here,levels(level).config);
        if run_solver
            generate_body_bmf(body_file,10,0,0,levels(level).N);
            timer = tic;
            results{level} = HydroMain(config_file);
            runtimes(level) = toc(timer);
        else
            cfg = read_config(config_file);
            if ~exist(cfg.files.results,'file')
                error('CRESTU:MissingConvergenceResult', ...
'Missing %s. Call Run_Convergence_Study(true).',cfg.files.results);
            end
            loaded = load(cfg.files.results,'results');
            results{level} = loaded.results;
        end
    end
    wamit_root = fullfile(root,'..','WAMIT');
    reference = struct();
    reference.IRR0 = load_reference(fullfile(wamit_root,'FullSphereIRR0'),results{1}.config);
    reference.IRR3 = load_reference(fullfile(wamit_root,'FullSphereIRR3'),results{1}.config);
    plot_file = fullfile(here,'Convergence_Study.png');
    drift_plot_file = fullfile(here,'MeanDrift_Validation.png');
    make_first_order_plot(results,levels,reference,plot_file);
    make_drift_plot(results{3},reference,drift_plot_file);
    summary = make_summary(results,levels,runtimes);
    study = struct('levels',levels,'results',{results},'runtimes',runtimes, ...
'reference',reference,'summary',summary,'plot_file',plot_file, ...
'drift_plot_file',drift_plot_file);
    writetable(summary,fullfile(here,'Convergence_Study.csv'));
    save(study_file,'study','-v7.3');
end

function reference = load_reference(folder,cfg)
% LOAD_REFERENCE Load first- and second-order WAMIT reference data.
%
% Syntax:
%   reference=load_reference(folder,cfg)
%
% Description:
%   Load first- and second-order WAMIT reference data.
%   The implementation preserves the CRESTU solver and data-field contracts.
%
% Inputs:
%   folder - Directory that contains reference files [-].
%   cfg - Validated CRESTU configuration structure [-].
%
% Outputs:
%   reference - WAMIT reference-data structure [-].
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

    reference.first_order = read_wamit_first_order(fullfile(folder,'hemisphere.1'),cfg.rho);
    reference.excitation = read_wamit_excitation(fullfile(folder,'hemisphere.2'),cfg.rho,cfg.grav,1,1);
    reference.rao = read_wamit_rao(fullfile(folder,'hemisphere.4'));
    reference.drift_momentum = read_wamit_mean_drift(fullfile(folder,'hemisphere.8'),cfg.rho,cfg.grav,1,1);
    reference.drift_pressure = read_wamit_mean_drift(fullfile(folder,'hemisphere.9'),cfg.rho,cfg.grav,1,1);
end

function make_first_order_plot(results,levels,reference,filename)
% MAKE_FIRST_ORDER_PLOT Plot first-order convergence and WAMIT comparisons.
%
% Syntax:
%   make_first_order_plot(results,levels,reference,filename)
%
% Description:
%   Plot first-order convergence and WAMIT comparisons.
%   The implementation preserves the CRESTU solver and data-field contracts.
%
% Inputs:
%   results - Hydrodynamic result structure or cell array [-].
%   levels - Grid-level definition structure [-].
%   reference - Reference-data structure [-].
%   filename - Output file path, character vector or string scalar [-].
%
% Outputs:
%   None.
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

    figure_handle = figure('Visible','off','Color','w','Position',[100,100,1200,820]);

%% Stage 2: Run the core calculation

    layout = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
    labels = {levels.name};
    colors = lines(3);
    axes_list = gobjects(4,1);
    for plotIndex = 1:4
        axes_list(plotIndex) = nexttile(layout);
        hold(axes_list(plotIndex),'on');
        grid(axes_list(plotIndex),'on');
    end
    for level = 1:3
        currentResult = results{level};
        omegaVector = currentResult.omegas;
        plot(axes_list(1),omegaVector,squeeze(currentResult.added_mass(3,3,:)),'o-','Color',colors(level,:),'DisplayName',labels{level});
        plot(axes_list(2),omegaVector,squeeze(currentResult.damping(3,3,:)),'o-','Color',colors(level,:),'DisplayName',labels{level});
        plot(axes_list(3),omegaVector,abs(squeeze(currentResult.excitation(3,1,:))),'o-','Color',colors(level,:),'DisplayName',labels{level});
        plot(axes_list(4),omegaVector,squeeze(currentResult.rao.amplitude(3,1,:)),'o-','Color',colors(level,:),'DisplayName',labels{level});
    end
    styles = {'k--','k:'};
    names = {'WAMIT IRR0','WAMIT IRR3'};
    refs = {reference.IRR0,reference.IRR3};
    for caseIndex = 1:2
        ref = refs{caseIndex};
        [wr,order] = sort(ref.first_order.omegas);
        plot(axes_list(1),wr,squeeze(ref.first_order.added_mass(3,3,order)),styles{caseIndex},'DisplayName',names{caseIndex});
        plot(axes_list(2),wr,squeeze(ref.first_order.damping(3,3,order)),styles{caseIndex},'DisplayName',names{caseIndex});
        [we,order] = sort(ref.excitation.omegas);
        plot(axes_list(3),we,abs(squeeze(ref.excitation.force(3,1,order))), ...
            styles{caseIndex},'DisplayName',names{caseIndex});
        [wx,order] = sort(ref.rao.omegas);
        plot(axes_list(4),wx,squeeze(ref.rao.amplitude(3,1,order)), ...
            styles{caseIndex},'DisplayName',names{caseIndex});
    end
    ylabel(axes_list(1),'A_{33} [kg]');
    ylabel(axes_list(2),'B_{33} [kg/s]');
    ylabel(axes_list(3),'|F_{z,exc}| [N/m]');
    ylabel(axes_list(4),'|\xi_3|/A');
    for plotIndex = 1:4
        xlabel(axes_list(plotIndex),'\omega [rad/s]');
        xlim(axes_list(plotIndex),[0.7,1.3]);
    end
    legend(axes_list(1),'Location','best');
    title(layout,'CRESTU hemisphere grid convergence');
    exportgraphics(figure_handle,filename,'Resolution',180);
    close(figure_handle);
end

function make_drift_plot(result,reference,filename)
% MAKE_DRIFT_PLOT Plot near-field and far-field mean-drift comparisons.
%
% Syntax:
%   make_drift_plot(result,reference,filename)
%
% Description:
%   Plot near-field and far-field mean-drift comparisons.
%   The implementation preserves the CRESTU solver and data-field contracts.
%
% Inputs:
%   result - Hydrodynamic result structure [-].
%   reference - Reference-data structure [-].
%   filename - Output file path, character vector or string scalar [-].
%
% Outputs:
%   None.
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

    if ~isfield(result,'drift') || ~result.drift.enabled
        return;
    end

%% Stage 2: Run the core calculation

    norm_scale = result.drift.normalization;
    omegaVector = result.omegas;
    near = squeeze(result.drift.near(1,1,:)) / norm_scale;
    far = squeeze(result.drift.far(1,1,:)) / norm_scale;
    fig = figure('Visible','off','Color','w');
    hold on;
    grid on;
    plot(omegaVector,near,'o-','DisplayName','CRESTU near field');
    plot(omegaVector,far,'s-','DisplayName','CRESTU far field');
    refs = {reference.IRR0,reference.IRR3};
    styles = {'k--','k:'};
    names = {'WAMIT IRR0 (.8/.9)','WAMIT IRR3 (.8/.9)'};
    for caseIndex = 1:2
        direct = refs{caseIndex}.drift_pressure;
        momentum = refs{caseIndex}.drift_momentum;
        [wd,order] = sort(direct.omegas);
        plot(wd,squeeze(direct.Cd(1,1,order)),styles{caseIndex}, ...
'DisplayName',[names{caseIndex},' pressure']);
        [wm,order] = sort(momentum.omegas);
        plot(wm,squeeze(momentum.Cd(1,1,order)),styles{caseIndex}, ...
'LineWidth',1.5,'DisplayName',[names{caseIndex},' momentum']);
    end
    xlabel('\omega [rad/s]');
    ylabel('C_{d,1}');
    xlim([0.7,1.3]);
    title('Mean surge drift: direct pressure vs Kochin momentum');
    legend('Location','best');
    exportgraphics(fig,filename,'Resolution',180);
    close(fig);
end

function summary = make_summary(results,levels,runtimes)
% MAKE_SUMMARY Assemble the convergence reconciliation table.
%
% Syntax:
%   summary=make_summary(results,levels,runtimes)
%
% Description:
%   Assemble the convergence reconciliation table.
%   The implementation preserves the CRESTU solver and data-field contracts.
%
% Inputs:
%   results - Hydrodynamic result structure or cell array [-].
%   levels - Grid-level definition structure [-].
%   runtimes - Elapsed solver times by grid level [s].
%
% Outputs:
%   summary - Reconciliation table for the completed study [-].
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

    rows = 0;

%% Stage 2: Run the core calculation

    for gridIndex = 1:3
        rows = rows + numel(results{gridIndex}.omegas);
    end
    grid_name = strings(rows,1);
    omega = zeros(rows,1);
    computational = zeros(rows,1);
    equivalent = zeros(rows,1);
    A33 = zeros(rows,1);
    B33 = zeros(rows,1);
    Fz = zeros(rows,1);
    heave_rao = zeros(rows,1);
    runtime = zeros(rows,1);
    cursor = 0;
    for caseIndex = 1:3
        currentResult = results{caseIndex};
        frequencyCount = numel(currentResult.omegas);
        idx = cursor + (1:frequencyCount);
        grid_name(idx) = levels(caseIndex).name;
        omega(idx) = currentResult.omegas;
        computational(idx) = currentResult.stats.total_dofs;
        equivalent(idx) = currentResult.symmetry.full_equivalent_dofs;
        A33(idx) = squeeze(currentResult.added_mass(3,3,:));
        B33(idx) = squeeze(currentResult.damping(3,3,:));
        Fz(idx) = abs(squeeze(currentResult.excitation(3,1,:)));
        heave_rao(idx) = squeeze(currentResult.rao.amplitude(3,1,:));
        runtime(idx) = runtimes(caseIndex);
        cursor = cursor + frequencyCount;
    end
    summary = table(grid_name,omega,computational,equivalent,A33,B33,Fz,heave_rao,runtime);
end
