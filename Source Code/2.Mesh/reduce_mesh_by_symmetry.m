function reduced = reduce_mesh_by_symmetry(mesh,isx,isy,tolerance)
% REDUCE_MESH_BY_SYMMETRY Execute the documented reduce_mesh_by_symmetry operation.
%
% Syntax:
%   reduced = reduce_mesh_by_symmetry(mesh,isx,isy,tolerance)
%
% Inputs:
%   mesh            : [struct] Boundary mesh with geometry expressed in SI units.
%   isx             : [logical scalar] Reflection-symmetry flag for the x = 0 plane.
%   isy             : [logical scalar] Reflection-symmetry flag for the y = 0 plane.
%   tolerance       : [documented value] Input required by the implemented function contract.
%
% Outputs:
%   reduced         : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
%REDUCE_MESH_BY_SYMMETRY Retain x>=0/y>=0 panels for geometric symmetry.
    if nargin<4||isempty(tolerance), tolerance=1e-9; end
    if ~isx&&~isy, reduced=mesh; return; end
    centers=reshape(mesh.centers,mesh.n_panels,3);
    vertices=reshape(mesh.vertices,mesh.n_panels,4,3);
    keep=true(mesh.n_panels,1);
    if isx, keep=keep&centers(:,1)>=-tolerance; end
    if isy, keep=keep&centers(:,2)>=-tolerance; end
    candidate=find(keep);
    if isx
        vx=reshape(vertices(candidate,:,1),numel(candidate),4);
        crossing=min(vx,[],2)<-tolerance&max(vx,[],2)>tolerance;
        if any(crossing)
            error('CRESTU:SymmetryPanelCrossing', ...
                '%d panels cross x=0; remesh with an x-symmetry edge.',nnz(crossing));
        end
    end
    if isy
        vy=reshape(vertices(candidate,:,2),numel(candidate),4);
        crossing=min(vy,[],2)<-tolerance&max(vy,[],2)>tolerance;
        if any(crossing)
            error('CRESTU:SymmetryPanelCrossing', ...
                '%d panels cross y=0; remesh with a y-symmetry edge.',nnz(crossing));
        end
    end
    fields={'vertices','centers','normals','areas','e1','e2','panel_type','mu_damping'};
    reduced=mesh;
    for k=1:numel(fields)
        name=fields{k};
        if isfield(mesh,name)&&~isempty(mesh.(name))
            value=mesh.(name); reduced.(name)=value(keep,:,:);
        end
    end
    reduced.n_panels=nnz(keep); reduced.isx=isx; reduced.isy=isy;
    reduced.header=sprintf('%s | reduced ISX=%d ISY=%d',mesh.header,isx,isy);
end
