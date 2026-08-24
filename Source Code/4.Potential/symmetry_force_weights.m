function weights = symmetry_force_weights(output_parity,solution_parity,isx,isy)
% SYMMETRY_FORCE_WEIGHTS Execute the documented symmetry_force_weights operation.
%
% Syntax:
%   weights = symmetry_force_weights(output_parity,solution_parity,isx,isy)
%
% Inputs:
%   output_parity   : [documented value] Input required by the implemented function contract.
%   solution_parity : [documented value] Input required by the implemented function contract.
%   isx             : [logical scalar] Reflection-symmetry flag for the x = 0 plane.
%   isy             : [logical scalar] Reflection-symmetry flag for the y = 0 plane.
%
% Outputs:
%   weights         : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
%SYMMETRY_FORCE_WEIGHTS Image-sum multiplier for reduced-surface integrals.
% output_parity: Nmode-by-2; solution_parity: Nsolution-by-2.
    weights=ones(size(output_parity,1),size(solution_parity,1));
    if isx, weights=weights.*(1+output_parity(:,1)*solution_parity(:,1).'); end
    if isy, weights=weights.*(1+output_parity(:,2)*solution_parity(:,2).'); end
end
