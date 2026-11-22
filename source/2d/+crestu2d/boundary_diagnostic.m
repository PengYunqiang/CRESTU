function [metric,detail] = boundary_diagnostic(geo,cfg,k,sigma,t,alpha)
%BOUNDARY_DIAGNOSTIC Check the solved density away from the collocation points.
% Report equal-point and one-point-per-panel arclength sampling estimates.
% The latter is NOT a quadrature-converged L2 norm; see round2 trace tests.
% Neither metric changes the solved coefficients or hides corner residuals.
    ids=cfg.dof_indices; scale=[1 1 cfg.normalization_length_m];D=diag(scale(ids));
    assert(~strcmp(cfg.solver.augmentation,'rigid_lid'),'CRESTU2D:TraceRepresentation', ...
        'This legacy trace helper accepts body and point-source strengths only; use saved solve diagnostics for rigid_lid. Omitting lid_strength would give a false trace.');
    if nargin<5,t=1/sqrt(3);end
    if nargin<6,alpha=[];end
    w=geo.c+t*geo.b;
    yz=[imag(w),cfg.waterline_z_m-real(w)];
    rhs=crestu2d.normal_modes(yz,geo.normal,cfg.reference_yz_m);rhs=rhs(:,ids);
    [~,sy,sd]=crestu2d.kernel(w,geo,k,cfg.solver.wave_order,cfg.fluid.domain,(1:geo.n)',t*ones(geo.n,1));
    actual=(geo.normal(:,1).*sy+geo.normal_d.*sd+pi*eye(geo.n))*sigma(:,1:numel(ids));
    if isfield(cfg.solver,'augmentation') && strcmp(cfg.solver.augmentation,'one_source_constraint')
        assert(~isempty(alpha),'CRESTU2D:Augmentation','Supply saved auxiliary strengths for an augmented trace.');
        [~,ay,ad]=crestu2d.point_green(k,imag(w),real(w),cfg.solver.auxiliary_source_y_m,0);
        actual=actual+(geo.normal(:,1).*ay+geo.normal_d.*ad)*alpha(1:numel(ids));
    end
    error=(actual-rhs)/D; target=rhs/D;
    err2=sum(abs(error).^2,2);rhs2=sum(abs(target).^2,2);
    metric.equal_point_relative=sqrt(sum(err2)/max(sum(rhs2),eps));
    metric.arclength_relative=sqrt(sum(geo.length.*err2)/max(sum(geo.length.*rhs2),eps));
    metric.maximum_scaled_absolute=max(sqrt(err2));
    metric.sample_local_t=t;
    metric.arclength_definition='one-point-per-panel estimate, not a converged L2 integral';
    metric.arclength_per_dof=sqrt(sum(geo.length.*abs(error).^2,1)./max(sum(geo.length.*abs(target).^2,1),eps));
    detail=table((1:geo.n)',yz(:,1),yz(:,2),geo.length,sqrt(err2),sqrt(rhs2),geo.length.*err2, ...
        'VariableNames',{'panel','y_m','z_m','length_m','residual_scaled','target_scaled','squared_error_arclength'});
    for j=1:numel(ids)
        detail.(['residual_' cfg.dofs{j}])=abs(error(:,j));
    end
end
