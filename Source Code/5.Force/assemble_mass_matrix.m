function mass_matrix = assemble_mass_matrix(mass_props)
% ASSEMBLE_MASS_MATRIX Execute the documented assemble_mass_matrix operation.
%
% Syntax:
%   mass_matrix = assemble_mass_matrix(mass_props)
%
% Inputs:
%   mass_props      : [struct array] Body mass, center of gravity, and inertia properties in SI units.
%
% Outputs:
%   mass_matrix     : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
%ASSEMBLE_MASS_MATRIX Build block-diagonal 6N rigid-body mass matrix at CGs.
    n=numel(mass_props); mass_matrix=zeros(6*n);
    for b=1:n
        idx=(b-1)*6+(1:6); m=mass_props(b).mass; inertia=mass_props(b).inertia;
        if ~isequal(size(inertia),[3,3])||m<=0||any(eig((inertia+inertia.')/2)<=0)
            error('CRESTU:MassProperties','Invalid mass properties for body %d.',b);
        end
        mass_matrix(idx,idx)=blkdiag(m*eye(3),(inertia+inertia.')/2);
    end
end
