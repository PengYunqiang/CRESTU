function result=hf_solve(nodes,opt)
%HF_SOLVE Unit-velocity HIGH-FREQUENCY radiation, not omega=0 or finite omega.
% Required: rho,depth,reference,quadrature_order. z=0 is the free surface.
    timer=tic;assert(opt.rho>0 && opt.depth>max(-nodes(:,2)));
    assert(size(nodes,2)==2 && all(isfinite(nodes),'all') && nodes(1,1)<nodes(end,1));
    assert(all(nodes([1 end],2)==0) && all(nodes(2:end-1,2)<0));
    delta=diff(nodes);len=vecnorm(delta,2,2);assert(all(len>0));
    signedArea=.5*sum(nodes(:,1).*circshift(nodes(:,2),-1)-nodes(:,2).*circshift(nodes(:,1),-1));assert(signedArea>0);
    nB=[delta(:,2),-delta(:,1)]./len;mid=(nodes(1:end-1,:)+nodes(2:end,:))/2;
    b=struct('mid',mid,'length',len,'tangent',delta./len,'normal',-nB,'n',numel(len));
    modes=crestu2d.normal_modes(mid,nB,opt.reference);
    [S,H]=crestu2d.hf_operators(mid,b,opt.depth,opt.quadrature_order);
    % <<<CORE>>> Exterior direct BIE: (pi I-H) phi + S qF=0, qF=-modes.
    M=pi*eye(b.n)-H;rhs=S*modes;
    rc=rcond(M);assert(rc>1e-13,'CRESTU2D:HFCondition','High-frequency direct system ill-conditioned.');
    phi=M\rhs;q=-modes;A=-opt.rho*(modes.'*(len.*phi));
    % <<< /CORE >>>
    L=nodes(end,1)-nodes(1,1);D=diag([1 1 L]);ah=D\A/D/(opt.rho*L^2);
    result=struct('schema','CRESTU-2D-HF-RESEARCH-1','method','DIRECT_BEM_EXACT_DIRICHLET_NEUMANN_STRIP_GREEN', ...
        'physical_limit','infinite_frequency','phase','exp(+i*omega*t)','normalization','unit velocity; roll unit angular velocity', ...
        'options',opt,'nodes_yz_m',nodes,'boundary',b,'body_normal',nB,'modes',modes,'potential',phi,'normal_derivative',q, ...
        'A',A,'A_hat',ah,'dofs',{{'sway','heave','roll'}},'units',{{'kg/m','kg/m','kg';'kg/m','kg/m','kg';'kg','kg','kg*m'}}, ...
        'linear_residual',norm(M*phi-rhs,'fro')/norm(rhs,'fro'),'matrix_rcond',rc, ...
        'reciprocity',norm(ah-ah.','fro')/norm(ah,'fro'),'minimum_scaled_eigenvalue',min(eig((ah+ah.')/2)), ...
        'elapsed_seconds',toc(timer),'scope','Global added mass only; no finite-omega pressure, damping, diffraction, or local-corner certification.');
    assert(all(isfinite(A),'all'),'CRESTU2D:HFNonfinite','Nonfinite high-frequency mass.');
end
