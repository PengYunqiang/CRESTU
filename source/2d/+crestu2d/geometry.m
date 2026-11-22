function geo = geometry(cfg)
%GEOMETRY Straight wetted panels, body outward normals in the transverse yz plane.
% Surface-piercing contour: left waterline -> keel -> right waterline.
% Closing it along the waterline makes a positive-area (CCW) body polygon.
    nodes=cfg.geometry.nodes_yz_m;
    assert(isnumeric(nodes) && isreal(nodes) && ismatrix(nodes) && size(nodes,2)==2 && size(nodes,1)>=4 && all(isfinite(nodes),'all'), ...
        'CRESTU2D:Geometry','Expected finite N-by-2 (y,z) nodes.');
    L=cfg.normalization_length_m; tol=1e-10*L;
    dz=nodes(:,2)-cfg.waterline_z_m;
    surfacePiercing=any(strcmp(cfg.fluid.domain,{'deep_water','finite_depth'}));
    if surfacePiercing
        assert(strcmp(cfg.geometry.orientation,'left_waterline_to_right_waterline'), ...
            'CRESTU2D:Orientation','Explicit left-to-right wetted contour required.');
        assert(nodes(1,1)<nodes(end,1),'CRESTU2D:Waterline', ...
            'Wetted contour must run left-to-right: reverse the node order explicitly; no automatic reversal is applied.');
        assert(max(abs(dz([1 end])))<tol,'CRESTU2D:Waterline', ...
            'Both endpoint elevations must equal waterline_z_m within the declared geometry tolerance.');
        assert(all(dz(2:end-1)<0),'CRESTU2D:Waterline', ...
            'All interior wetted nodes must be strictly submerged; do not add a waterline closure panel.');
    else
        assert(strcmp(cfg.geometry.orientation,'counterclockwise_closed'), ...
            'CRESTU2D:Orientation','Unbounded-fluid orientation must be counterclockwise_closed.');
        assert(norm(nodes(1,:)-nodes(end,:))<tol,'CRESTU2D:Orientation', ...
            'Unbounded-fluid contour must be explicitly closed by repeating its first node.');
    end
    if strcmp(cfg.fluid.domain,'finite_depth')
        assert(all(dz>-cfg.fluid.depth_m),'CRESTU2D:BottomClearance','Every body node must lie strictly above the seabed.');
        assert(cfg.solver.side_extent_m>max(abs(nodes(:,1))),'CRESTU2D:SideBoundary','Side boundaries must lie outside the entire section.');
    end
    a=nodes(1:end-1,:); b=nodes(2:end,:); delta=b-a;
    len=sqrt(sum(delta.^2,2));
    assert(all(len>tol),'CRESTU2D:Geometry','Repeated node or zero-length panel.');
    area=.5*sum(nodes(:,1).*circshift(nodes(:,2),-1)-nodes(:,2).*circshift(nodes(:,1),-1));
    assert(area>tol^2,'CRESTU2D:Orientation','Negative/zero body area: reverse the node order explicitly.');
    if cfg.geometry.check_self_intersection
        closed=[nodes;nodes(1,:)];
        if norm(closed(end,:)-closed(end-1,:))<tol, closed(end,:)=[]; end
        n=size(closed,1)-1;
        for i=1:n
            for j=i+2:n
                if i==1 && j==n, continue; end
                p=closed(i,:); r=closed(i+1,:)-p;
                q=closed(j,:); s=closed(j+1,:)-q;
                det=cross2(r,s);
                if abs(det)>eps*L^2
                    t=cross2(q-p,s)/det; u=cross2(q-p,r)/det;
                    assert(~(t>=0 && t<=1 && u>=0 && u<=1),'CRESTU2D:Geometry','Self-intersecting contour.');
                elseif abs(cross2(q-p,r))<tol*norm(r)
                    t1=dot(q-p,r)/dot(r,r);t2=dot(q+s-p,r)/dot(r,r);
                    assert(max(min(t1,t2),0)>min(max(t1,t2),1),'CRESTU2D:Geometry','Overlapping contour segments.');
                end
            end
        end
    end
    geo.nodes=nodes; geo.mid=(a+b)/2; geo.length=len;
    geo.normal=[delta(:,2),-delta(:,1)]./len;
    geo.area_m2=area; geo.n=numel(len);
    % w=d+i*y; y origin is fixed and only z is shifted to the waterline.
    w=cfg.waterline_z_m-nodes(:,2)+1i*nodes(:,1);
    geo.c=(w(1:end-1)+w(2:end))/2; geo.b=(w(2:end)-w(1:end-1))/2;
    geo.jac=len/2;
    geo.normal_d=-geo.normal(:,2);
    geo.modes=crestu2d.normal_modes(geo.mid,geo.normal,cfg.reference_yz_m);
    % Diagnostics describe the supplied straight-panel geometry only. A turn
    % in a polygonal approximation does not establish an analytic corner.
    tangent=delta./len;
    if surfacePiercing
        vertex=(2:geo.n)'; incoming=tangent(1:end-1,:); outgoing=tangent(2:end,:);
    else
        vertex=(1:geo.n)'; incoming=circshift(tangent,1,1); outgoing=tangent;
    end
    turn=atan2(incoming(:,1).*outgoing(:,2)-incoming(:,2).*outgoing(:,1),sum(incoming.*outgoing,2));
    geo.diagnostics=struct('orientation',cfg.geometry.orientation,'signed_area_m2',area, ...
        'minimum_panel_length_m',min(len),'maximum_panel_length_m',max(len), ...
        'self_intersection_check_enabled',logical(cfg.geometry.check_self_intersection), ...
        'vertex_indices',vertex,'vertex_yz_m',nodes(vertex,:), ...
        'vertex_turn_angle_rad',turn,'fluid_side_vertex_angle_rad',pi+turn, ...
        'corner_interpretation','Angles belong to the input polygon; no smooth-curve or physical-corner classification is inferred.');
    if surfacePiercing
        wetTangent=[tangent(1,:);-tangent(end,:)];
        geo.diagnostics.waterline_endpoint_error_m=dz([1 end]);
        geo.diagnostics.waterline_y_m=nodes([1 end],1);
        geo.diagnostics.endpoint_into_wet_tangent_yz=wetTangent;
        geo.diagnostics.endpoint_into_wet_angle_rad=atan2(wetTangent(:,2),wetTangent(:,1));
    end
end
function z=cross2(x,y),z=x(1)*y(2)-x(2)*y(1);end
