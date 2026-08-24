function [A,S_body,assembly_info] = assemble_rankine_matrix(N,centers,normals,vertices,stats,omega,g,cfg,parity)
% ASSEMBLE_RANKINE_MATRIX Execute the documented assemble_rankine_matrix operation.
%
% Syntax:
%   [A,S_body,assembly_info] = assemble_rankine_matrix(N,centers,normals,vertices,stats,omega,g,cfg,parity)
%
% Inputs:
%   N               : [integer scalar or array] Discrete count or index required by the algorithm.
%   centers         : [N x 3] Panel collocation points, in m.
%   normals         : [N x 3] Unit panel normals pointing into the fluid.
%   vertices        : [N x 4 x 3] Quadrilateral panel vertices, in m.
%   stats           : [documented value] Input required by the implemented function contract.
%   omega           : [scalar] Angular frequency, in rad/s.
%   g               : [scalar] Gravitational acceleration, in m/s^2.
%   cfg             : [struct] Validated CRESTU configuration with SI-valued physical fields.
%   parity          : [K x 2] Reflection parity signs for x and y symmetry planes.
%
% Outputs:
%   A               : [documented value] Function result; dimensions and units follow the implemented contract.
%   S_body          : [documented value] Function result; dimensions and units follow the implemented contract.
%   assembly_info   : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
%ASSEMBLE_RANKINE_MATRIX Mixed-boundary Rankine BIE with symmetry images.
% parity=[px,py], where +1 is symmetric and -1 antisymmetric.  Image source
% contributions enforce Neumann/Dirichlet conditions on omitted symmetry
% planes, so no panels are required on x=0 or y=0.
    if nargin<9||isempty(parity), parity=[1,1]; end
    validateattributes(N,{'numeric'},{'scalar','integer','positive'});
    validateattributes(omega,{'numeric'},{'scalar','real','positive','finite'});
    validateattributes(g,{'numeric'},{'scalar','real','positive','finite'});
    parity=reshape(parity,1,2);
    if any(abs(parity)~=1), error('CRESTU:SymmetryParity','Parity entries must be +1 or -1.'); end
    centers=reshape(centers,N,3); normals=reshape(normals,N,3); vertices=reshape(vertices,N,4,3);
    n_body=stats.total_body_panels; n_fs=stats.fs_panels; n_sb=stats.seabed_panels;
    n_ff=N-n_body-n_fs-n_sb;
    if n_ff<0||(isfield(stats,'farfield_panels')&&stats.farfield_panels~=n_ff)
        error('CRESTU:DomainCounts','Panel statistics do not match N=%d.',N);
    end
    if any(~isfinite(centers(:)))||any(~isfinite(normals(:)))||any(~isfinite(vertices(:)))
        error('CRESTU:NonfiniteGeometry','BEM geometry contains NaN or Inf.');
    end
    [image_vertices,image_normals,image_weights,image_labels]=prepare_images( ...
        vertices,normals,cfg.isx,cfg.isy,parity);

    nu=zeros(N,1); [k0,~]=solve_wave_dispersion(omega,g,cfg.water_depth);
    sponge_width=cfg.fs.r_outer-cfg.fs.r_inner;
    effective_mu0=min(2.5,max(cfg.fs.mu0,12/(k0*sponge_width)));
    fs_idx=n_body+(1:n_fs);
    if n_fs>0
        radii=sqrt(sum(centers(fs_idx,1:2).^2,2)); damping=zeros(n_fs,1);
        mask=radii>cfg.fs.r_inner;
        damping(mask)=effective_mu0*((radii(mask)-cfg.fs.r_inner)/sponge_width).^2;
        damping=max(0,min(effective_mu0,damping));
        % exp(i*w*t): exp(-i*k*r), hence Im(k)<0 gives spatial decay.
        nu(fs_idx)=k0*(1-1i*damping);
    end
    ff_idx=n_body+n_fs+n_sb+(1:n_ff); nu(ff_idx)=k0;
    A=complex(zeros(N,N)); S_body=complex(zeros(N,n_body)); inv4pi=1/(4*pi);
    fprintf('>>> Rankine assembly N=%d, parity=[%+d,%+d], images=%d\n', ...
        N,parity(1),parity(2),numel(image_weights)); timer=tic;
    for i=1:N
        rowA=complex(zeros(1,N)); rowS=complex(zeros(1,n_body)); x=centers(i,:);
        for j=1:N
            Gsum=0; Dsum=0;
            for im=1:numel(image_weights)
                panel=reshape(image_vertices(j,:,:,im),4,3);
                normal=reshape(image_normals(j,:,im),1,3);
                [G,dGdn]=rankine_panel_integrals(panel,x,normal);
                Gsum=Gsum+image_weights(im)*G;
                Dsum=Dsum+image_weights(im)*dGdn;
            end
            Sij=Gsum*inv4pi; Dij=-Dsum*inv4pi;
            if i==j, Dij=Dij+0.5; end
            if j<=n_body
                rowA(j)=Dij; rowS(j)=Sij;
            elseif j<=n_body+n_fs
                rowA(j)=Dij-nu(j)*Sij;
            elseif j<=n_body+n_fs+n_sb
                rowA(j)=Dij;
            else
                % Far-field normals point inward: dphi/dn=+i*k*phi.
                rowA(j)=Dij-1i*nu(j)*Sij;
            end
        end
        A(i,:)=rowA; S_body(i,:)=rowS;
    end
    elapsed=toc(timer);
    assembly_info=struct('elapsed_seconds',elapsed,'parity',parity, ...
        'image_count',numel(image_weights),'image_labels',{image_labels}, ...
        'wavenumber',k0,'wavelength',2*pi/k0,'effective_mu0',effective_mu0, ...
        'sponge_width',sponge_width,'n_unknowns',N);
    fprintf('>>> Rankine assembly completed in %.3f s.\n',elapsed);
end

function [all_vertices,all_normals,weights,labels]=prepare_images(vertices,normals,isx,isy,parity)
% PREPARE_IMAGES Execute the documented prepare_images operation.
%
% Syntax:
%   [all_vertices,all_normals,weights,labels]=prepare_images(vertices,normals,isx,isy,parity)
%
% Inputs:
%   vertices        : [N x 4 x 3] Quadrilateral panel vertices, in m.
%   normals         : [N x 3] Unit panel normals pointing into the fluid.
%   isx             : [logical scalar] Reflection-symmetry flag for the x = 0 plane.
%   isy             : [logical scalar] Reflection-symmetry flag for the y = 0 plane.
%   parity          : [K x 2] Reflection parity signs for x and y symmetry planes.
%
% Outputs:
%   all_vertices    : [documented value] Function result; dimensions and units follow the implemented contract.
%   all_normals     : [documented value] Function result; dimensions and units follow the implemented contract.
%   weights         : [documented value] Function result; dimensions and units follow the implemented contract.
%   labels          : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
    flags=[0,0];
    if isx, flags=[flags;1,0]; end
    if isy, flags=[flags;0,1]; end
    if isx&&isy, flags=[flags;1,1]; end
    n=size(vertices,1); ni=size(flags,1);
    all_vertices=zeros(n,4,3,ni); all_normals=zeros(n,3,ni);
    weights=ones(ni,1); labels=cell(ni,1);
    for q=1:ni
        v=vertices; normal=normals;
        if flags(q,1), v(:,:,1)=-v(:,:,1); normal(:,1)=-normal(:,1); end
        if flags(q,2), v(:,:,2)=-v(:,:,2); normal(:,2)=-normal(:,2); end
        if mod(sum(flags(q,:)),2)==1, v(:,[2,4],:)=v(:,[4,2],:); end
        all_vertices(:,:,:,q)=v; all_normals(:,:,q)=normal;
        weights(q)=parity(1)^flags(q,1)*parity(2)^flags(q,2);
        labels{q}=sprintf('reflect_x=%d,reflect_y=%d',flags(q,1),flags(q,2));
    end
end
