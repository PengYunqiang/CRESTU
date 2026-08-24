function nj = compute_generalized_normals(centers, normals, body_list)
% COMPUTE_GENERALIZED_NORMALS Execute the documented compute_generalized_normals operation.
%
% Syntax:
%   nj = compute_generalized_normals(centers, normals, body_list)
%
% Inputs:
%   centers         : [N x 3] Panel collocation points, in m.
%   normals         : [N x 3] Unit panel normals pointing into the fluid.
%   body_list       : [documented value] Input required by the implemented function contract.
%
% Outputs:
%   nj              : [N x 6Nb] Generalized panel normals for all body modes.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
%COMPUTE_GENERALIZED_NORMALS Build the panel-by-6N generalized-normal map.
% Translational columns are n. Rotational columns are (r-XG) x n.
    centers=normalize_xyz(centers,'centers'); normals=normalize_xyz(normals,'normals');
    if size(centers,1)~=size(normals,1)
        error('CRESTU:GeometrySize','Centers and normals must have the same row count.');
    end
    if isstruct(body_list), body_list=num2cell(body_list); end
    if ~iscell(body_list)||isempty(body_list)
        error('CRESTU:BodyList','body_list must be a nonempty cell/struct array.');
    end
    n_bodies=numel(body_list); counts=zeros(n_bodies,1);
    for b=1:n_bodies
        if ~isfield(body_list{b},'n_panels')||~isfield(body_list{b},'cg')
            error('CRESTU:BodyMetadata','Body %d requires n_panels and cg fields.',b);
        end
        counts(b)=body_list{b}.n_panels;
    end
    if sum(counts)~=size(centers,1)
        error('CRESTU:BodyPanelCount','Body panel counts total %d, geometry has %d rows.', ...
            sum(counts),size(centers,1));
    end
    nj=zeros(size(centers,1),6*n_bodies); first=1;
    for b=1:n_bodies
        idx=first:(first+counts(b)-1); cols=(b-1)*6+(1:6);
        n_body=normals(idx,:); cg=reshape(body_list{b}.cg,1,3);
        nj(idx,cols(1:3))=n_body;
        nj(idx,cols(4:6))=cross(centers(idx,:)-cg,n_body,2);
        first=idx(end)+1;
    end
end

function xyz=normalize_xyz(xyz,label)
% NORMALIZE_XYZ Execute the documented normalize_xyz operation.
%
% Syntax:
%   xyz=normalize_xyz(xyz,label)
%
% Inputs:
%   xyz             : [documented value] Input required by the implemented function contract.
%   label           : [documented value] Input required by the implemented function contract.
%
% Outputs:
%   xyz             : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
    if isempty(xyz), xyz=zeros(0,3); end
    if isvector(xyz)&&numel(xyz)==3, xyz=reshape(xyz,1,3); end
    if ~isnumeric(xyz)||~ismatrix(xyz)||size(xyz,2)~=3||any(~isfinite(xyz(:)))
        error('CRESTU:GeometryShape','%s must be a finite N-by-3 numeric array.',label);
    end
end
