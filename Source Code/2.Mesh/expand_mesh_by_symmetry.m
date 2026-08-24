function full_mesh = expand_mesh_by_symmetry(mesh,isx,isy)
% EXPAND_MESH_BY_SYMMETRY Execute the documented expand_mesh_by_symmetry operation.
%
% Syntax:
%   full_mesh = expand_mesh_by_symmetry(mesh,isx,isy)
%
% Inputs:
%   mesh            : [struct] Boundary mesh with geometry expressed in SI units.
%   isx             : [logical scalar] Reflection-symmetry flag for the x = 0 plane.
%   isy             : [logical scalar] Reflection-symmetry flag for the y = 0 plane.
%
% Outputs:
%   full_mesh       : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
%EXPAND_MESH_BY_SYMMETRY Reflect a reduced mesh for auxiliary mesh creation.
% This utility is geometric only; potential parity signs are handled by the
% image kernel and are not stored in this mesh.
    full_mesh=mesh;
    if isx, full_mesh=append_reflection(full_mesh,1); end
    if isy, full_mesh=append_reflection(full_mesh,2); end
    full_mesh.isx=0; full_mesh.isy=0;
    full_mesh.header=sprintf('%s | geometric full expansion',mesh.header);
end

function out=append_reflection(in,axis_index)
% APPEND_REFLECTION Execute the documented append_reflection operation.
%
% Syntax:
%   out=append_reflection(in,axis_index)
%
% Inputs:
%   in              : [documented value] Input required by the implemented function contract.
%   axis_index      : [documented value] Input required by the implemented function contract.
%
% Outputs:
%   out             : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
    reflected=in; reflected.vertices(:,:,axis_index)=-reflected.vertices(:,:,axis_index);
    reflected.centers(:,axis_index)=-reflected.centers(:,axis_index);
    reflected.normals(:,axis_index)=-reflected.normals(:,axis_index);
    reflected.e1(:,axis_index)=-reflected.e1(:,axis_index);
    reflected.e2(:,axis_index)=-reflected.e2(:,axis_index);
    % Reflection reverses vertex orientation. Swap 2 and 4 to preserve the
    % original normal represented by the reflected normal array.
    reflected.vertices(:,[2,4],:)=reflected.vertices(:,[4,2],:);
    out=in; fields={'vertices','centers','normals','areas','e1','e2','panel_type','mu_damping'};
    for k=1:numel(fields)
        name=fields{k};
        if isfield(in,name)&&~isempty(in.(name))
            out.(name)=cat(1,in.(name),reflected.(name));
        end
    end
    out.n_panels=size(out.centers,1);
end
