function mesh_global = transform_body_mesh(mesh_local, body_cfg)
% TRANSFORM_BODY_MESH Execute the documented transform_body_mesh operation.
%
% Syntax:
%   mesh_global = transform_body_mesh(mesh_local, body_cfg)
%
% Inputs:
%   mesh_local      : [struct] Body mesh in its local coordinate frame.
%   body_cfg        : [documented value] Input required by the implemented function contract.
%
% Outputs:
%   mesh_global     : [struct] Body mesh transformed to the global coordinate frame.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
%

    mesh_global = mesh_local;
    x0 = body_cfg.pos;
    yaw_rad = body_cfg.yaw * pi / 180.0;

    Rz = [ cos(yaw_rad), -sin(yaw_rad),  0;
           sin(yaw_rad),  cos(yaw_rad),  0;
           0,             0,             1 ];

    np = mesh_local.n_panels;

    for v = 1:4
        pts_loc = reshape(mesh_local.vertices(:, v, :), np, 3);
        pts_glob = pts_loc * Rz' + repmat(x0, [np, 1]);
        mesh_global.vertices(:, v, :) = pts_glob;
    end

    mesh_global.centers = mesh_local.centers * Rz' + repmat(x0, [np, 1]);
    mesh_global.normals = mesh_local.normals * Rz';
    mesh_global.e1      = mesh_local.e1 * Rz';
    mesh_global.e2      = mesh_local.e2 * Rz';

    if isfield(mesh_local, 'hydrostatics') && mesh_local.hydrostatics.V_mean > 0
        cb_loc = mesh_local.hydrostatics.center_of_buoyancy;
        mesh_global.hydrostatics.center_of_buoyancy = cb_loc * Rz' + x0;
    end
end
