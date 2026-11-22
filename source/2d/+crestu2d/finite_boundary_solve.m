function out=finite_boundary_solve(nodes,omega,opt)
%FINITE_BOUNDARY_SOLVE Direct Laplace BEM with bottom images and analytic DtN.
% Required opt: rho,g,depth,side,reference,normalization_length,
% surface_step,bottom_step,side_panels,direction_y. Unit elevation amplitude.
    timer=tic;rho=opt.rho;g=opt.g;depth=opt.depth;side=opt.side;
    assert(strcmp(opt.bottom_representation,'image_neumann'), ...
        'CRESTU2D:FiniteBottom','The admitted finite solver uses the exact Neumann image kernel.');
    direction=opt.direction_y;assert(ismember(direction,[-1 1]));
    ref=opt.reference(:).';L=opt.normalization_length;nu=omega^2/g;
    assert(isreal(nodes) && size(nodes,2)==2 && all(isfinite(nodes),'all') && size(nodes,1)>=4);
    assert(omega>0 && rho>0 && g>0 && depth>max(-nodes(:,2)) && side>max(abs(nodes(:,1))));
    assert(all(nodes([1 end],2)==0) && all(nodes(2:end-1,2)<0) && nodes(1,1)<nodes(end,1));
    assert(opt.surface_step>0 && opt.bottom_step>0 && opt.side_panels>=2 && opt.side_panels==fix(opt.side_panels));
    b=make_boundary(nodes,opt);
    [S,H]=crestu2d.finite_operators(b.mid,b);
    constantError=norm((pi*eye(b.n)-H)*ones(b.n,1),inf);
    assert(constantError<1e-9,'CRESTU2D:FiniteOrientation','Closed-domain constant identity failed; check BIE signs and normals.');
    Q=zeros(b.n);Q(b.fs,b.fs)=nu*eye(numel(b.fs));
    right=crestu2d.finite_modes(omega,g,depth,b.mid(b.right,2));
    left=crestu2d.finite_modes(omega,g,depth,b.mid(b.left,2));
    Q(b.right,b.right)=right.T;Q(b.left,b.left)=left.T;
    body=b.body;nb=-b.normal(body,:);xy=b.mid(body,:);
    modes=[nb,(xy(:,1)-ref(1)).*nb(:,2)-(xy(:,2)-ref(2)).*nb(:,1)];
    k=right.k;[shape,shapez]=incident(k,xy(:,2),depth);
    phii=1i*g/omega*shape.*exp(-1i*k*direction*xy(:,1));
    gradIy=-1i*k*direction*phii;gradIz=1i*g/omega*shapez.*exp(-1i*k*direction*xy(:,1));
    imposed=complex(zeros(b.n,4));imposed(body,1:3)=-modes;
    imposed(body,4)=nb(:,1).*gradIy+nb(:,2).*gradIz;
    % <<<CORE>>> Direct Green identity, all normals outward from fluid.
    matrix=pi*eye(b.n)-H+S*Q;rhs=-S*imposed;
    rc=rcond(matrix);assert(rc>1e-13,'CRESTU2D:FiniteSystem','Direct finite-domain system is ill-conditioned.');
    phi=matrix\rhs;q=Q*phi+imposed;
    C=modes.'*(b.length(body).*phi(body,1:3));
    A=-rho*real(C);B=rho*omega*imag(C);
    F=1i*rho*omega*(modes.'*(b.length(body).*(phi(body,4)+phii)));
    FK=1i*rho*omega*(modes.'*(b.length(body).*phii));
    % <<< /CORE >>>
    rightCoefficients=right.V\phi(b.right,:);leftCoefficients=left.V\phi(b.left,:);
    far=[rightCoefficients(1,1:3);leftCoefficients(1,1:3)];
    Benergy=rho*omega*k*right.propagating_mode_norm_m*real(far'*far);
    globalFar=far*exp(1i*k*side);
    Fhaskind=2i*rho*g*k*right.propagating_mode_norm_m*globalFar(1+(direction==1),:).';
    D=diag([1 1 L]);
    out=struct('method','FINITE_DEPTH_DIRECT_IMAGE_BEM_ANALYTIC_MODAL_DTN', ...
        'phase','exp(+i*omega*t)','omega',omega,'options',opt,'nodes',nodes, ...
        'A',A,'B',B,'F',F,'F_FK',FK,'A_hat',D\A/D/(rho*L^2), ...
        'B_hat',D\B/D/(rho*omega*L^2),'F_hat',D\F/(rho*g*L), ...
        'boundary',b,'potential',phi,'normal_derivative',q, ...
        'body_radiation_potential',phi(body,1:3),'body_scattered_potential_per_amplitude',phi(body,4), ...
        'body_incident_potential_per_amplitude',phii,'body_modes',modes, ...
        'right_modes',right,'left_modes',left,'farfield_radiation_potential',far, ...
        'farfield_radiation_potential_global',globalFar,'F_Haskind',Fhaskind, ...
        'haskind_relative',norm(D\(F-Fhaskind))/max(norm(D\Fhaskind),eps), ...
        'equation_residual',norm(matrix*phi-rhs,'fro')/max(norm(rhs,'fro'),eps), ...
        'constant_identity_error',constantError,'matrix_rcond',rc, ...
        'reciprocity_A',norm(D\(A-A.')/D,'fro')/max(norm(D\A/D,'fro'),eps), ...
        'reciprocity_B',norm(D\(B-B.')/D,'fro')/max(norm(D\B/D,'fro'),eps), ...
        'energy_relative',norm(D\(B-Benergy)/D,'fro')/max(norm(D\Benergy/D,'fro'),eps), ...
        'B_energy',Benergy,'elapsed_seconds',toc(timer));
    % Reconstruct a genuine fluid-side gradient, independent of q=0 imposed
    % in the algebraic bottom rows. Do not add a boundary jump to this probe.
    bp=b.bottom_probe_geometry;shift=1e-4*bp.length;
    probes=bp.mid+(bp.length/(2*sqrt(3))).*bp.tangent-shift.*bp.normal;
    [~,gy,gz]=crestu2d.finite_field(out,probes);
    residual=bp.normal(:,1).*gy+bp.normal(:,2).*gz;
    velocityScale=[1 1 L g*k/omega];scaled=residual./velocityScale;
    out.bottom_probe=struct('points_yz_m',probes,'normal_velocity',residual, ...
        'rms_scaled',sqrt(sum(bp.length.*sum(abs(scaled).^2,2))/sum(bp.length)), ...
        'offset_fraction',1e-4,'normalization',velocityScale);
    assert(all(isfinite([A(:);B(:);F(:);phi(:)])),'CRESTU2D:FiniteNonfinite','Nonfinite finite-depth BEM result.');
    fprintf('Finite direct BEM N=%d body=%d side_modes=%d H=%.5g k=%.6g residual=%.3g bottom=%.3g\n', ...
        b.n,numel(body),opt.side_panels,depth,k,out.equation_residual,out.bottom_probe.rms_scaled);
end
function b=make_boundary(nodes,opt)
    side=opt.side;depth=opt.depth;
    nleft=ceil((nodes(1,1)+side)/opt.surface_step);
    nright=ceil((side-nodes(end,1))/opt.surface_step);
    nz=opt.side_panels;nbottom=ceil(2*side/opt.bottom_step);
    pieces={ [linspace(-side,nodes(1,1),nleft+1).',zeros(nleft+1,1)],nodes, ...
        [linspace(nodes(end,1),side,nright+1).',zeros(nright+1,1)], ...
        [side*ones(nz+1,1),linspace(0,-depth,nz+1).'], ...
        [linspace(side,-side,nbottom+1).',-depth*ones(nbottom+1,1)], ...
        [-side*ones(nz+1,1),linspace(-depth,0,nz+1).']};
    a=[];c=[];groups=cell(1,6);
    for j=1:6
        if j==5,groups{j}=[];continue;end
        p=pieces{j};groups{j}=size(a,1)+(1:size(p,1)-1);
        a=[a;p(1:end-1,:)];c=[c;p(2:end,:)]; %#ok<AGROW>
    end
    delta=c-a;len=sqrt(sum(delta.^2,2));assert(all(len>0));
    bottom=pieces{5};db=diff(bottom);lb=sqrt(sum(db.^2,2));
    bottomProbe=struct('mid',(bottom(1:end-1,:)+bottom(2:end,:))/2,'length',lb, ...
        'tangent',db./lb,'normal',[-db(:,2),db(:,1)]./lb);
    b=struct('a',a,'b',c,'mid',(a+c)/2,'length',len,'tangent',delta./len, ...
        'normal',[-delta(:,2),delta(:,1)]./len,'n',numel(len), ...
        'fs',[groups{1} groups{3}],'body',groups{2},'right',groups{4}, ...
        'bottom',groups{5},'left',groups{6},'orientation','clockwise; fluid on right; outward normal on left', ...
        'bottom_probe_geometry',bottomProbe,'depth',depth,'bottom_representation',opt.bottom_representation);
end
function [f,fz]=incident(k,z,depth)
    expterm=exp(-2*k*(z+depth));den=1+exp(-2*k*depth);
    f=exp(k*z).*(1+expterm)/den;fz=k*exp(k*z).*(1-expterm)/den;
end
