function [C,details] = compute_hydrostatic_matrix(body_list,cfg)
% COMPUTE_HYDROSTATIC_MATRIX Execute the documented compute_hydrostatic_matrix operation.
%
% Syntax:
%   [C,details] = compute_hydrostatic_matrix(body_list,cfg)
%
% Inputs:
%   body_list       : [documented value] Input required by the implemented function contract.
%   cfg             : [struct] Validated CRESTU configuration with SI-valued physical fields.
%
% Outputs:
%   C               : [documented value] Function result; dimensions and units follow the implemented contract.
%   details         : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
%COMPUTE_HYDROSTATIC_MATRIX Assemble block-diagonal 6N restoring matrix.
% Waterplane moments are evaluated about each body's CG horizontal axes.
    if isstruct(body_list), body_list=num2cell(body_list); end
    n=numel(body_list); C=zeros(6*n);
    template=struct('body_id',0,'waterplane_area',0,'center_of_flotation',zeros(1,2), ...
        'Ixx_wp',0,'Iyy_wp',0,'Ixy_wp',0,'first_moment_x',0,'first_moment_y',0, ...
        'displaced_volume',0,'center_of_buoyancy',zeros(1,3),'block',zeros(6));
    details=repmat(template,n,1);
    for b=1:n
        mesh=body_list{b}; cg=reshape(cfg.mass_props(b).cg,1,3);
        wl=extract_waterline(mesh,cfg.z_tol); xy=wl.nodes;
        if cfg.isx||cfg.isy
            wl=complete_waterline_by_symmetry(wl,cfg.isx,cfg.isy); xy=wl.nodes;
        end
        [Awp,centroid,Ixx,Iyy,Ixy,Qx,Qy]=polygon_properties(xy,cg(1:2));
        hs=mesh.hydrostatics;
        volume=hs.V_mean;
        if ~(isfinite(volume)&&volume>0), volume=hs.Vz; end
        if ~(isfinite(volume)&&volume>0)
            error('CRESTU:HydrostaticVolume','Body %d has no positive displaced volume.',b);
        end
        zb=hs.center_of_buoyancy(3); mass=cfg.mass_props(b).mass; block=zeros(6);
        rg=cfg.rho*cfg.grav;
        block(3,3)=rg*Awp;
        block(3,4)=rg*Qy; block(4,3)=block(3,4);
        block(3,5)=-rg*Qx; block(5,3)=block(3,5);
        vertical_term=rg*volume*zb-mass*cfg.grav*cg(3);
        block(4,4)=rg*Ixx+vertical_term;
        block(5,5)=rg*Iyy+vertical_term;
        block(4,5)=-rg*Ixy; block(5,4)=block(4,5);
        idx=(b-1)*6+(1:6); C(idx,idx)=block;
        details(b)=struct('body_id',cfg.bodies(b).id,'waterplane_area',Awp, ...
            'center_of_flotation',centroid,'Ixx_wp',Ixx,'Iyy_wp',Iyy,'Ixy_wp',Ixy, ...
            'first_moment_x',Qx,'first_moment_y',Qy,'displaced_volume',volume, ...
            'center_of_buoyancy',hs.center_of_buoyancy,'block',block);
    end
end

function [A,centroid,Ixx,Iyy,Ixy,Qx,Qy]=polygon_properties(xy,reference)
% POLYGON_PROPERTIES Execute the documented polygon_properties operation.
%
% Syntax:
%   [A,centroid,Ixx,Iyy,Ixy,Qx,Qy]=polygon_properties(xy,reference)
%
% Inputs:
%   xy              : [documented value] Input required by the implemented function contract.
%   reference       : [documented value] Input required by the implemented function contract.
%
% Outputs:
%   A               : [documented value] Function result; dimensions and units follow the implemented contract.
%   centroid        : [documented value] Function result; dimensions and units follow the implemented contract.
%   Ixx             : [documented value] Function result; dimensions and units follow the implemented contract.
%   Iyy             : [documented value] Function result; dimensions and units follow the implemented contract.
%   Ixy             : [documented value] Function result; dimensions and units follow the implemented contract.
%   Qx              : [documented value] Function result; dimensions and units follow the implemented contract.
%   Qy              : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
    if size(xy,1)<3, error('CRESTU:WaterplanePolygon','Waterline needs at least three vertices.'); end
    xy=xy-reference; x=xy(:,1); y=xy(:,2); xn=x([2:end,1]); yn=y([2:end,1]);
    cross_terms=x.*yn-xn.*y; signed_A=0.5*sum(cross_terms);
    if abs(signed_A)<eps, error('CRESTU:WaterplaneArea','Waterplane polygon area is zero.'); end
    s=sign(signed_A); A=abs(signed_A);
    Qx=s*sum((x+xn).*cross_terms)/6; Qy=s*sum((y+yn).*cross_terms)/6;
    Ixx=s*sum((y.^2+y.*yn+yn.^2).*cross_terms)/12;
    Iyy=s*sum((x.^2+x.*xn+xn.^2).*cross_terms)/12;
    Ixy=s*sum((2*x.*y+x.*yn+xn.*y+2*xn.*yn).*cross_terms)/24;
    centroid=reference+[Qx/A,Qy/A];
end
