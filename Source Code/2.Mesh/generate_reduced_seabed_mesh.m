function mesh = generate_reduced_seabed_mesh(cfg,full_circumference_nodes)
% GENERATE_REDUCED_SEABED_MESH Execute the documented generate_reduced_seabed_mesh operation.
%
% Syntax:
%   mesh = generate_reduced_seabed_mesh(cfg,full_circumference_nodes)
%
% Inputs:
%   cfg             : [struct] Validated CRESTU configuration with SI-valued physical fields.
%   full_circumference_nodes: [documented value] Input required by the implemented function contract.
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
%GENERATE_REDUCED_SEABED_MESH Polar half/quarter disk with symmetry edges.
    if cfg.water_depth<=0, error('CRESTU:SeabedDepth','Finite positive depth required.'); end
    if cfg.isx&&cfg.isy
        theta0=0; theta1=pi/2; divisor=4;
    elseif cfg.isx, theta0=-pi/2; theta1=pi/2; divisor=2;
    elseif cfg.isy, theta0=0; theta1=pi; divisor=2;
    else, theta0=0; theta1=2*pi; divisor=1;
    end
    ntheta=max(4,round(full_circumference_nodes/divisor));
    theta=linspace(theta0,theta1,ntheta+1);
    nr=max(4,cfg.fs.nr_near+cfg.fs.nr_sponge); radial=cfg.fs.r_outer*(linspace(0,1,nr+1).^0.8);
    np=nr*ntheta; vertices=zeros(np,4,3); p=0; z=-cfg.water_depth;
    for i=1:nr
        for j=1:ntheta
            p=p+1; r0=radial(i); r1=radial(i+1); t0=theta(j); t1=theta(j+1);
            vertices(p,1,:)=[r0*cos(t0),r0*sin(t0),z];
            vertices(p,2,:)=[r1*cos(t0),r1*sin(t0),z];
            vertices(p,3,:)=[r1*cos(t1),r1*sin(t1),z];
            vertices(p,4,:)=[r0*cos(t1),r0*sin(t1),z];
        end
    end
    [centers,areas,e1,e2]=panel_geometry(vertices);
    mesh=struct('header',sprintf('Reduced polar seabed ISX=%d ISY=%d',cfg.isx,cfg.isy), ...
        'ulen',1,'panel_type',repmat(4,np,1),'isx',cfg.isx,'isy',cfg.isy, ...
        'n_panels',np,'vertices',vertices,'centers',centers,'normals',repmat([0,0,1],np,1), ...
        'areas',areas,'e1',e1,'e2',e2,'hydrostatics',empty_hydro());
end

function [centers,areas,e1,e2]=panel_geometry(vertices)
% PANEL_GEOMETRY Execute the documented panel_geometry operation.
%
% Syntax:
%   [centers,areas,e1,e2]=panel_geometry(vertices)
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
        P=reshape(vertices(p,:,:),4,3); centers(p,:)=mean(P,1);
        areas(p)=0.5*norm(cross(P(3,:)-P(1,:),P(4,:)-P(2,:)));
        tangent=P(2,:)-P(1,:);
        if norm(tangent)<1e-12, tangent=P(3,:)-P(1,:); end
        e1(p,:)=tangent/norm(tangent); e2(p,:)=cross([0,0,1],e1(p,:));
    end
end
function h=empty_hydro()
% EMPTY_HYDRO Execute the documented empty_hydro operation.
%
% Syntax:
%   h=empty_hydro()
%
% Inputs:
%   None.
%
% Outputs:
%   h               : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
    h=struct('Vx',0,'Vy',0,'Vz',0,'V_mean',0,'center_of_buoyancy',[0,0,0]);
end
