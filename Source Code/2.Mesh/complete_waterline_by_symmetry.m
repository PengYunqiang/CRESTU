function full_waterline = complete_waterline_by_symmetry(waterline,isx,isy)
% COMPLETE_WATERLINE_BY_SYMMETRY Execute the documented complete_waterline_by_symmetry operation.
%
% Syntax:
%   full_waterline = complete_waterline_by_symmetry(waterline,isx,isy)
%
% Inputs:
%   waterline       : [struct] Ordered waterline nodes in the z = 0 plane, in m.
%   isx             : [logical scalar] Reflection-symmetry flag for the x = 0 plane.
%   isy             : [logical scalar] Reflection-symmetry flag for the y = 0 plane.
%
% Outputs:
%   full_waterline  : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
%COMPLETE_WATERLINE_BY_SYMMETRY Reconstruct a full contour from half/quarter.
    points=waterline.nodes;
    if isx, points=[points;[-points(:,1),points(:,2)]]; end
    if isy, points=[points;[points(:,1),-points(:,2)]]; end
    points=uniquetol(points,1e-9,'ByRows',true);
    center=mean(points,1); angles=atan2(points(:,2)-center(2),points(:,1)-center(1));
    [~,order]=sort(angles); points=points(order,:);
    full_waterline=struct('nodes',points,'theta',atan2(points(:,2),points(:,1)), ...
        'r',sqrt(sum(points.^2,2)),'n_pts',size(points,1),'is_closed',true);
end
