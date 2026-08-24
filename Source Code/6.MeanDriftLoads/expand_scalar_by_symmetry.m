function expanded = expand_scalar_by_symmetry(values,parity,isx,isy)
% EXPAND_SCALAR_BY_SYMMETRY Execute the documented expand_scalar_by_symmetry operation.
%
% Syntax:
%   expanded = expand_scalar_by_symmetry(values,parity,isx,isy)
%
% Inputs:
%   values          : [documented value] Input required by the implemented function contract.
%   parity          : [K x 2] Reflection parity signs for x and y symmetry planes.
%   isx             : [logical scalar] Reflection-symmetry flag for the x = 0 plane.
%   isy             : [logical scalar] Reflection-symmetry flag for the y = 0 plane.
%
% Outputs:
%   expanded        : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
%EXPAND_SCALAR_BY_SYMMETRY Reconstruct scalar fields on reflected panels.
% values is Npanel-by-Ncomponent; parity is Ncomponent-by-[x,y].
    if isempty(values), expanded=values; return; end
    parity=reshape(parity,size(values,2),2);
    flags=[0,0];
    if isx, flags=[flags;1,0]; end
    if isy, flags=[flags;0,1]; end
    if isx&&isy, flags=[0,0;1,0;0,1;1,1]; end
    expanded=complex(zeros(size(values,1)*size(flags,1),size(values,2)));
    for q=1:size(flags,1)
        rows=(q-1)*size(values,1)+(1:size(values,1));
        weights=(parity(:,1).^flags(q,1)).*(parity(:,2).^flags(q,2));
        expanded(rows,:)=values.*weights.';
    end
end
