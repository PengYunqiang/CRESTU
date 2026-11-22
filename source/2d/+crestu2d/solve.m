function result = solve(cfg)
%SOLVE Independent 2-D source BEM. No 3-D kernels, defaults or meshes used.
    if strcmp(cfg.fluid.domain,'finite_depth'),result=crestu2d.solve_finite(cfg);return;end
    geo=crestu2d.geometry(cfg); ids=cfg.dof_indices;
    rho=cfg.fluid.rho_kg_m3; g=cfg.fluid.g_m_s2; L=cfg.normalization_length_m;
    nd=numel(ids); nf=numel(cfg.omega_rad_s); scale=[1 1 L]; D=diag(scale(ids));
    diffraction=strcmp(cfg.problem,'radiation_diffraction');
    augmented=strcmp(cfg.solver.augmentation,'one_source_constraint');
    withLid=strcmp(cfg.solver.augmentation,'rigid_lid');constrained=augmented || withLid;
    cap=[];
    if withLid
        ends=linspace(geo.nodes(1,1),geo.nodes(end,1),cfg.solver.lid_panels+1)';
        cap=struct('c',1i*(ends(1:end-1)+ends(2:end))/2,'b',1i*diff(ends)/2, ...
            'jac',diff(ends)/2,'n',cfg.solver.lid_panels,'endpoints_y_m',ends);
    end
    if augmented
        sourceY=cfg.solver.auxiliary_source_y_m;p=cfg.solver.interior_constraint_yz_m;
        assert(sourceY>geo.nodes(1,1) && sourceY<geo.nodes(end,1), ...
            'CRESTU2D:Augmentation','Auxiliary source must lie strictly inside the internal waterline segment.');
        [inside,on]=inpolygon(p(1),p(2),geo.nodes(:,1),geo.nodes(:,2));
        assert(inside && ~on && p(2)<cfg.waterline_z_m, ...
            'CRESTU2D:Augmentation','Constraint point must lie strictly inside the closed body polygon.');
        interiorW=cfg.waterline_z_m-p(2)+1i*p(1);
    end
    [t,wt]=crestu2d.gauss(cfg.solver.force_order);
    wf=geo.c+geo.b*t.'; wf=wf(:);
    weights=geo.jac*wt.'; weights=weights(:);
    own=repmat((1:geo.n)',numel(t),1); local=kron(t,ones(geo.n,1));
    normal=repmat(geo.normal,numel(t),1);
    yz=[imag(wf),cfg.waterline_z_m-real(wf)];
    modes=crestu2d.normal_modes(yz,normal,cfg.reference_yz_m); modes=modes(:,ids);
    [tq,wq]=crestu2d.gauss(cfg.solver.wave_order); ws=geo.c+geo.b*tq.';
    % Separate off-collocation boundary check at a fixed interior point.
    tc=1/sqrt(3); wc=geo.c+tc*geo.b;
    checkmode=crestu2d.normal_modes([imag(wc),cfg.waterline_z_m-real(wc)],geo.normal,cfg.reference_yz_m);
    result=struct('schema','CRESTU-2D-RESULT-1','version','0.3-research', ...
        'name',cfg.name,'input',cfg,'dofs',{cfg.dofs},'omega_rad_s',cfg.omega_rad_s, ...
        'reference_yz_m',cfg.reference_yz_m,'geometry',geo, ...
        'phase','exp(+i*omega*t)','radiation_normalization','unit velocity; roll unit angular velocity', ...
        'capability','Deep-water linear 2D; explicit rigid lid suppresses the continuous interior spectrum. Discrete qualification is bounded by saved geometry/frequency/grid evidence.');
    result.augmentation=cfg.solver.augmentation;
    result.auxiliary_strength=complex(zeros(double(augmented),nd+diffraction,nf));
    result.internal_lid=cap;
    if withLid,result.lid_strength=complex(zeros(cap.n,nd+diffraction,nf));end
    result.A=zeros(nd,nd,nf); result.B=result.A;
    result.A_hat=result.A; result.B_hat=result.A;
    result.F_exc_per_amplitude=complex(zeros(nd,diffraction*nf));
    result.F_FK_per_amplitude=result.F_exc_per_amplitude;
    result.farfield_velocity_potential=complex(zeros(2,nd,nf));
    result.radiation_potential=complex(zeros(numel(wf),nd,nf));
    result.source_strength=complex(zeros(geo.n,nd+diffraction,nf));
    result.force_points_yz_m=yz; result.force_weights_m=weights;
    result.units=struct('A',{{'kg/m','kg/m','kg';'kg/m','kg/m','kg';'kg','kg','kg*m'}}, ...
        'B',{{'kg/(m*s)','kg/(m*s)','kg/s';'kg/(m*s)','kg/(m*s)','kg/s';'kg/s','kg/s','kg*m/s'}}, ...
        'F_exc',{{'N/m','N/m','N'}},'F_exc_per_amplitude',{{'N/m^2','N/m^2','N/m'}}, ...
        'note','Units tables use full sway/heave/roll order; all loads per unit longitudinal length.');
    result.normalization='A_hat=D^-1*A*D^-1/(rho*L^2); B_hat=D^-1*B*D^-1/(rho*omega*L^2); D=diag(1,1,L) selected DOFs';
    for f=1:nf
        omega=cfg.omega_rad_s(f); k=omega^2/g;
        [~,sy,sd]=crestu2d.kernel(geo.c,geo,k,cfg.solver.wave_order,cfg.fluid.domain,(1:geo.n)',zeros(geo.n,1));
        % <<<CORE>>> Exterior normal-derivative jump is +pi for G~+log(r).
        M=geo.normal(:,1).*sy+geo.normal_d.*sd+pi*eye(geo.n);
        rhs=geo.modes(:,ids);
        if diffraction
            inc=1i*exp(-k*real(geo.c)-1i*k*cfg.incident.direction_y*imag(geo.c));
            dninc=inc.*(-1i*k*cfg.incident.direction_y*geo.normal(:,1)-k*geo.normal_d);
            rhs=[rhs,-dninc]; %#ok<AGROW>
        end
        baseRcond=rcond(M);system=M;systemRhs=rhs;
        if augmented
            [~,ay,ad]=crestu2d.point_green(k,imag(geo.c),real(geo.c),sourceY,0);
            gn=geo.normal(:,1).*ay+geo.normal_d.*ad;
            gi=crestu2d.kernel(interiorW,geo,k,cfg.solver.wave_order,'deep_water');
            ai=crestu2d.point_green(k,imag(interiorW),real(interiorW),sourceY,0);
            % The extra zero is a representation constraint, not a fluid BC.
            % Solve for beta=alpha/L; scale the potential row by L as well.
            system=[M L*gn;gi/L ai];systemRhs=[rhs;zeros(1,nd+diffraction)];
        elseif withLid
            [~,ly,ld]=crestu2d.kernel(geo.c,cap,k,cfg.solver.wave_order,'deep_water');
            sib=crestu2d.kernel(cap.c,geo,k,cfg.solver.wave_order,'deep_water');
            sii=crestu2d.lid_potential(imag(cap.c),cap.endpoints_y_m,k,cfg.solver.lid_quadrature_order);
            % Interior rigid lid: d_z(phi)=k*phi-2*pi*lambda=0.
            system=[M geo.normal(:,1).*ly+geo.normal_d.*ld;-k*sib 2*pi*eye(cap.n)-k*sii];
            systemRhs=[rhs;zeros(cap.n,nd+diffraction)];
        end
        if constrained
            rs=vecnorm(system,2,2);mr=system./rs;cs=vecnorm(mr,2,1);scaled=mr./cs;
            rc=rcond(scaled);
        else
            rc=baseRcond;
        end
        assert(rc>1e-13,'CRESTU2D:Singular','Ill-conditioned source system (rcond %.3g); possible irregular frequency.',rc);
        if constrained
            strength=(scaled\(systemRhs./rs))./cs.';
            sigma=strength(1:geo.n,:);
            if augmented
                alpha=L*strength(end,:);result.auxiliary_strength(:,:,f)=alpha;
            else
                lambda=strength(geo.n+1:end,:);result.lid_strength(:,:,f)=lambda;
            end
        else
            strength=M\rhs;sigma=strength;
        end
        sf=crestu2d.kernel(wf,geo,k,cfg.solver.wave_order,cfg.fluid.domain,own,local);
        phi=sf*sigma;
        if augmented,phi=phi+crestu2d.point_green(k,imag(wf),real(wf),sourceY,0)*alpha;end
        if withLid,phi=phi+crestu2d.kernel(wf,cap,k,cfg.solver.wave_order,'deep_water')*lambda;end
        rad=phi(:,1:nd);
        I=modes.'*(weights.*rad);
        A=-rho*real(I); B=rho*omega*imag(I);
        % No symmetrization, clipping, artificial damping or coefficient fitting.
        result.A(:,:,f)=A; result.B(:,:,f)=B;
        ah=(D\A/D)/(rho*L^2); bh=(D\B/D)/(rho*omega*L^2);
        result.A_hat(:,:,f)=ah; result.B_hat(:,:,f)=bh;
        result.radiation_potential(:,:,f)=rad; result.source_strength(:,:,f)=sigma;
        metric=struct('omega_rad_s',omega,'k_m_inv',k,'matrix_rcond',rc, ...
            'unaugmented_matrix_rcond',baseRcond,'augmentation_active',constrained, ...
            'linear_residual',norm(system*strength-systemRhs,'fro')/max(norm(systemRhs,'fro'),eps), ...
            'reciprocity_A',norm(ah-ah.','fro')/max(norm(ah,'fro'),eps), ...
            'reciprocity_B',norm(bh-bh.','fro')/max(norm(bh,'fro'),eps), ...
            'min_damping_eigenvalue_hat',min(eig((bh+bh.')/2)));
        [~,cy,cd]=crestu2d.kernel(wc,geo,k,cfg.solver.wave_order,cfg.fluid.domain,(1:geo.n)',tc*ones(geo.n,1));
        cv=(geo.normal(:,1).*cy+geo.normal_d.*cd+pi*eye(geo.n))*sigma(:,1:nd);
        if augmented
            [~,ay,ad]=crestu2d.point_green(k,imag(wc),real(wc),sourceY,0);
            cv=cv+(geo.normal(:,1).*ay+geo.normal_d.*ad)*alpha(1:nd);
        elseif withLid
            [~,ly,ld]=crestu2d.kernel(wc,cap,k,cfg.solver.wave_order,'deep_water');
            cv=cv+(geo.normal(:,1).*ly+geo.normal_d.*ld)*lambda(:,1:nd);
        end
        metric.offcollocation_bc_relative=norm((cv-checkmode(:,ids))/D,'fro')/max(norm(checkmode(:,ids)/D,'fro'),eps);
        if strcmp(cfg.fluid.domain,'deep_water')
            right=(exp(-k*real(ws)+1i*k*imag(ws))*wq).*geo.jac;
            left=(exp(-k*real(ws)-1i*k*imag(ws))*wq).*geo.jac;
            h=2i*pi*[right.';left.']*sigma(:,1:nd);
            if augmented,h=h+2i*pi*exp(1i*k*[sourceY;-sourceY])*alpha(1:nd);end
            if withLid
                wl=cap.c+cap.b*tq.';
                hr=(exp(-k*conj(wl))*wq).*cap.jac;hl=(exp(-k*wl)*wq).*cap.jac;
                h=h+2i*pi*[hr.';hl.']*lambda(:,1:nd);
            end
            result.farfield_velocity_potential(:,:,f)=h;
            benergy=rho*omega/2*real(h'*h);
            metric.energy_relative=norm(D\(B-benergy)/D,'fro')/max(norm(D\benergy/D,'fro'),eps);
            metric.farfield_hermitian_imag_relative=norm(imag(h'*h),'fro')/max(norm(h'*h,'fro'),eps);
            result.B_energy(:,:,f)=benergy;
        end
        if diffraction
            dir=cfg.incident.direction_y;
            uinc=1i*exp(-k*real(wf)-1i*k*dir*imag(wf));
            total=uinc+phi(:,end);
            F=1i*rho*g*(modes.'*(weights.*total));
            Ffk=1i*rho*g*(modes.'*(weights.*uinc));
            dninc=uinc.*(-1i*k*dir*normal(:,1)+k*normal(:,2));
            Frecip=1i*rho*g*(modes.'*(weights.*uinc)-rad.'*(weights.*dninc));
            incoming_side=1+(dir==1); Fhaskind=1i*rho*g*h(incoming_side,:).';
            result.F_exc_per_amplitude(:,f)=F; result.F_FK_per_amplitude(:,f)=Ffk;
            result.F_Haskind_per_amplitude(:,f)=Fhaskind;
            result.diffraction_potential(:,f)=phi(:,end);
            metric.excitation_reciprocity_relative=norm(D\(F-Frecip))/max(norm(D\Frecip),eps);
            metric.haskind_relative=norm(D\(F-Fhaskind))/max(norm(D\Fhaskind),eps);
        end
        % <<< /CORE >>>
        assert(all(isfinite([A(:);B(:);phi(:)])),'CRESTU2D:Nonfinite','Nonfinite numerical result.');
        % Diagnostic review flags are not automatic physical certification.
        metric.boundary_review_required=metric.offcollocation_bc_relative>.05;
        metric.conditioning_review_required=rc<1e-3;
        result.diagnostics(f)=metric;
        fprintf('CRESTU-2D %s: omega=%.8g rad/s, N=%d, rcond=%.3g, residual=%.3g\n', ...
            cfg.name,omega,geo.n,rc,metric.linear_residual);
    end
    result.F_exc=result.F_exc_per_amplitude;
    if diffraction, result.F_exc=result.F_exc*cfg.incident.amplitude_m; end
end
