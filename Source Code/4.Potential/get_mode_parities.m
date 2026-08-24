function parity = get_mode_parities(n_bodies,isx,isy)
% GET_MODE_PARITIES Execute the documented get_mode_parities operation.
%
% Syntax:
%   parity = get_mode_parities(n_bodies,isx,isy)
%
% Inputs:
%   n_bodies        : [integer scalar or array] Discrete count or index required by the algorithm.
%   isx             : [logical scalar] Reflection-symmetry flag for the x = 0 plane.
%   isy             : [logical scalar] Reflection-symmetry flag for the y = 0 plane.
%
% Outputs:
%   parity          : [K x 2] Reflection parity signs for x and y symmetry planes.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
%GET_MODE_PARITIES Return [x-plane,y-plane] parity for all 6N modes.
% +1 is symmetric (Neumann on plane); -1 is antisymmetric (Dirichlet).
    px=[-1,1,1,1,-1,-1]; py=[1,-1,1,-1,1,-1];
    parity=[repmat(px(:),n_bodies,1),repmat(py(:),n_bodies,1)];
    if ~isx, parity(:,1)=1; end
    if ~isy, parity(:,2)=1; end
end
