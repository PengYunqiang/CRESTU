function report = audit_rankine_wave_resolution(domain, omega)
% AUDIT_RANKINE_WAVE_RESOLUTION Inspect ACTUAL generated boundary spacing.
% This is a screening diagnostic, not a convergence certificate or a gate.
% Radial and total-edge sizes are reported separately. A large azimuthal
% edge in an axisymmetric far field need not imply a large radial phase error.
    cfg = domain.cfg;
    [k,~] = solve_wave_dispersion(omega,cfg.grav,cfg.water_depth,cfg.numerics);
    lambda = 2*pi/k;
    report = struct('schemaVersion',1,'omegaRadPerS',omega, ...
        'wavenumberPerM',k,'wavelengthM',lambda, ...
        'screeningOnly',true,'convergenceCertified',false);
    names = {'fs','seabed','farfield'};
    for j=1:numel(names)
        name=names{j}; m=domain.(name);
        if isempty(m)
            report.(name)=struct('panelCount',0); continue
        end
        v=m.vertices; nv=v(:,[2,3,4,1],:);
        edges=sqrt(sum((nv-v).^2,3));
        radius=sqrt(v(:,:,1).^2+v(:,:,2).^2);
        radial=abs(radius(:,[2,3,4,1])-radius);
        report.(name)=struct('panelCount',m.n_panels, ...
            'maxEdgeM',max(edges(:)),'medianEdgeM',median(edges(:)), ...
            'maxRadialEdgeM',max(radial(:)), ...
            'maxRadialEdgeOverLambda',max(radial(:))/lambda, ...
            'maxEdgeOverLambda',max(edges(:))/lambda, ...
            'maxSqrtAreaM',max(sqrt(m.areas)));
    end
    report.bodyPanelCount=domain.stats.total_body_panels;
    report.configuredOuterRadiusM=cfg.fs.r_outer;
    outerVertices=domain.fs.vertices;
    if ~isempty(domain.farfield),outerVertices=domain.farfield.vertices;end
    report.outerRadiusM=max(hypot(outerVertices(:,:,1),outerVertices(:,:,2)),[],'all');
    report.spongeStartRadiusM=cfg.fs.r_inner;
    report.finiteDepthVerticalFactor=NaN;
    if cfg.water_depth>0
        kh=k*cfg.water_depth;
        report.finiteDepthVerticalFactor=2*exp(-kh)/(1+exp(-2*kh));
    end
    if report.fs.panelCount>0 && report.fs.maxRadialEdgeOverLambda>1/8
        warning('CRESTU:UnresolvedFreeSurface', ...
            ['Actual FS max radial edge is %.3g lambda at omega=%.5g. ', ...
             'Refine the FS and test convergence; a small matrix residual is insufficient.'], ...
             report.fs.maxRadialEdgeOverLambda,omega);
    end
    if cfg.water_depth>0 && report.seabed.panelCount>0 && ...
            report.finiteDepthVerticalFactor>0.1 && ...
            report.seabed.maxRadialEdgeOverLambda>1/8
        warning('CRESTU:UnresolvedFiniteDepthBottom', ...
            ['Actual bottom max radial edge is %.3g lambda; ', ...
             'the propagating-wave bottom factor is %.3g. Refine the bottom independently.'], ...
            report.seabed.maxRadialEdgeOverLambda,report.finiteDepthVerticalFactor);
    end
end
