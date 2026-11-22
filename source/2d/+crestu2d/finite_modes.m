function modes=finite_modes(omega,g,depth,z)
%FINITE_MODES Analytic vertical eigenmodes, one per side trace DOF.
% exp(+i omega t): outward DtN eigenvalues are -i*k and -kappa_n.
    nu=omega^2/g;count=numel(z);z=z(:);
    assert(omega>0 && g>0 && depth>0 && count>=2 && all(z>=-depth & z<=0));
    hi=max(nu+1/depth,2*sqrt(nu/depth));
    while hi*tanh(hi*depth)<nu,hi=2*hi;end
    k=fzero(@(q)q*tanh(q*depth)-nu,[0 hi]);
    kappas=zeros(count-1,1);
    for j=1:count-1
        lo=(j-.5)*pi;up=j*pi;offset=16*eps(up);
        x=fzero(@(x)x*tan(x)+nu*depth,[lo+offset up-offset]);
        kappas(j)=x/depth;
    end
    shape=exp(k*z).*(1+exp(-2*k*(z+depth)))/(1+exp(-2*k*depth));
    V=[shape,cos((z+depth)*kappas.')];
    lambda=[-1i*k;-kappas];
    assert(rcond(V)>1e-12,'CRESTU2D:FiniteModes','Sampled analytic modal basis is singular.');
    T=(V.*lambda.')/V;
    fs=[k*tanh(k*depth)-nu; -kappas.*sin(kappas*depth)-nu*cos(kappas*depth)];
    bottomZ=-depth;
    bottom=[k*exp(k*bottomZ)*(1-exp(-2*k*(bottomZ+depth)))/(1+exp(-2*k*depth)); ...
        -kappas.*sin(kappas*(bottomZ+depth))];
    cg=omega/(2*k)*(1+2*k*depth/sinh(2*k*depth));
    norm0=depth/2*(1-tanh(k*depth)^2)+tanh(k*depth)/(2*k);
    modes=struct('k',k,'kappa',kappas,'V',V,'T',T,'outward_eigenvalues',lambda, ...
        'number_of_modes',count,'group_velocity',cg,'propagating_mode_norm_m',norm0,'basis_rcond',rcond(V), ...
        'flux_identity_relative',abs((omega/g*cg)/(k*norm0)-1), ...
        'dispersion_relative',abs(k*tanh(k*depth)/nu-1), ...
        'surface_mode_residual',max(abs(fs)./(nu+[k;kappas])), ...
        'bottom_mode_derivative',max(abs(bottom)), ...
        'dtn_mode_residual',norm(T*V-V.*lambda.','fro')/norm(V.*lambda.','fro'), ...
        'scope','Complete square sampled analytic basis: one propagating and Ns-1 evanescent modes; Ns remains a convergence parameter.');
end
