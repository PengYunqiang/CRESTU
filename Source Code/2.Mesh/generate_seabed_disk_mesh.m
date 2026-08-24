function mesh = generate_seabed_disk_mesh(cfg,mesh_fs)
% GENERATE_SEABED_DISK_MESH Execute the documented generate_seabed_disk_mesh operation.
%
% Syntax:
%   mesh = generate_seabed_disk_mesh(cfg,mesh_fs)
%
% Inputs:
%   cfg             : [struct] Validated CRESTU configuration with SI-valued physical fields.
%   mesh_fs         : [documented value] Input required by the implemented function contract.
%
% Outputs:
%   mesh            : [struct] Boundary mesh with geometry expressed in SI units.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
%GENERATE_SEABED_DISK_MESH Close a multi-body domain with a full flat disk.
    v=reshape(mesh_fs.vertices,mesh_fs.n_panels,4,3);
    r=sqrt(v(:,:,1).^2+v(:,:,2).^2); mask=abs(r-cfg.fs.r_outer)<1e-6*max(1,cfg.fs.r_outer);
    xy=[reshape(v(:,:,1),[],1),reshape(v(:,:,2),[],1)]; boundary=xy(mask(:),:);
    if size(boundary,1)<8
        theta=linspace(0,2*pi,65).'; theta(end)=[]; boundary=cfg.fs.r_outer*[cos(theta),sin(theta)];
    else
        angles=atan2(boundary(:,2),boundary(:,1)); [~,order]=sort(angles);
        boundary=unique(boundary(order,:),'rows','stable');
    end
    n=size(boundary,1); constraints=[(1:n).',[2:n,1].']; spacing=max(cfg.fs.r_outer/14,0.75);
    a=-cfg.fs.r_outer:spacing:cfg.fs.r_outer; [gx,gy]=meshgrid(a,a); pts=[gx(:),gy(:)];
    pts=pts(sum(pts.^2,2)<(0.985*cfg.fs.r_outer)^2,:); points=[boundary;pts];
    [points,~,map]=uniquetol(points,1e-10,'ByRows',true); constraints=map(constraints);
    constraints=constraints(constraints(:,1)~=constraints(:,2),:);
    dt=delaunayTriangulation(points,constraints); tri=dt.ConnectivityList; points=dt.Points;
    c=(points(tri(:,1),:)+points(tri(:,2),:)+points(tri(:,3),:))/3;
    tri=tri(sum(c.^2,2)<=cfg.fs.r_outer^2,:); np=size(tri,1); vertices=zeros(np,4,3);
    for q=1:3, vertices(:,q,1:2)=points(tri(:,q),:); end
    vertices(:,4,:)=vertices(:,3,:); vertices(:,:,3)=-cfg.water_depth;
    [centers,areas,e1,e2]=triangle_geometry(vertices);
    mesh=struct('header',sprintf('Multi-body disk seabed h=%.3g',cfg.water_depth),'ulen',1, ...
        'panel_type',repmat(4,np,1),'isx',0,'isy',0,'n_panels',np,'vertices',vertices, ...
        'centers',centers,'normals',repmat([0,0,1],np,1),'areas',areas,'e1',e1,'e2',e2, ...
        'hydrostatics',struct('Vx',0,'Vy',0,'Vz',0,'V_mean',0,'center_of_buoyancy',[0,0,0]));
end

function [centers,areas,e1,e2]=triangle_geometry(vertices)
% TRIANGLE_GEOMETRY Execute the documented triangle_geometry operation.
%
% Syntax:
%   [centers,areas,e1,e2]=triangle_geometry(vertices)
%
% Inputs:
%   vertices        : [N x 4 x 3] Quadrilateral panel vertices, in m.
%
% Outputs:
%   centers         : [N x 3] Panel collocation points, in m.
%   areas           : [N x 1] Panel areas, in m^2.
%   e1              : [documented value] Function result; dimensions and units follow the implemented contract.
%   e2              : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
    np=size(vertices,1); centers=zeros(np,3); areas=zeros(np,1); e1=zeros(np,3); e2=zeros(np,3);
    for p=1:np
        p1=reshape(vertices(p,1,:),1,3); p2=reshape(vertices(p,2,:),1,3); p3=reshape(vertices(p,3,:),1,3);
        centers(p,:)=(p1+p2+p3)/3; areas(p)=0.5*norm(cross(p2-p1,p3-p1));
        e1(p,:)=(p2-p1)/norm(p2-p1); e2(p,:)=cross([0,0,1],e1(p,:));
    end
end
