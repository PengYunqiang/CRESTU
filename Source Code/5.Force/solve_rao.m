function rao = solve_rao(omegas,added_mass,damping,hydrostatic,excitation,cfg)
% SOLVE_RAO Execute the documented solve_rao operation.
%
% Syntax:
%   rao = solve_rao(omegas,added_mass,damping,hydrostatic,excitation,cfg)
%
% Inputs:
%   omegas          : [1 x Nf] Angular frequencies, in rad/s.
%   added_mass      : [Ndof x Ndof x Nf] Added-mass matrices in SI translational/rotational units.
%   damping         : [Ndof x Ndof x Nf] Radiation-damping matrices in SI translational/rotational units.
%   hydrostatic     : [Ndof x Ndof] Hydrostatic restoring matrix in SI units.
%   excitation      : [Ndof x Nh x Nf] Complex first-order excitation loads in SI units.
%   cfg             : [struct] Validated CRESTU configuration with SI-valued physical fields.
%
% Outputs:
%   rao             : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
%SOLVE_RAO Solve [-w^2(M+A)+i*w*B+C] Xi=Fexc for all wave cases.
% Arrays A/B are ndof-by-ndof-by-nfreq; excitation is ndof-by-nhead-by-nfreq.
    omegas=reshape(omegas,1,[]); nf=numel(omegas); ndof=6*cfg.n_bodies;
    added_mass=normalize_coeff(added_mass,ndof,nf,'added mass');
    damping=normalize_coeff(damping,ndof,nf,'damping');
    if ~isequal(size(hydrostatic),[ndof,ndof])
        error('CRESTU:HydrostaticShape','Hydrostatic matrix must be %d-by-%d.',ndof,ndof);
    end
    if ismatrix(excitation) && nf==1, excitation=reshape(excitation,ndof,size(excitation,2),1); end
    if size(excitation,1)~=ndof||size(excitation,3)~=nf
        error('CRESTU:ExcitationShape','Excitation must be ndof-by-nheading-by-nfrequency.');
    end
    nh=size(excitation,2); M=assemble_mass_matrix(cfg.mass_props);
    xi=complex(zeros(ndof,nh,nf)); dynamic_matrix=complex(zeros(ndof,ndof,nf)); rconds=zeros(nf,1);
    for k=1:nf
        w=omegas(k); Z=-w^2*(M+added_mass(:,:,k))+1i*w*damping(:,:,k)+hydrostatic;
        dynamic_matrix(:,:,k)=Z; rconds(k)=rcond(Z);
        if rconds(k)<1e-12, warning('CRESTU:IllConditionedRAO','RAO matrix at omega=%g has rcond=%g.',w,rconds(k)); end
        xi(:,:,k)=Z\excitation(:,:,k);
    end
    rao=struct('complex',xi,'amplitude',abs(xi),'phase_deg',angle(xi)*180/pi, ...
        'mass_matrix',M,'dynamic_matrix',dynamic_matrix,'rcond',rconds, ...
        'omegas',omegas,'headings',cfg.wave.headings);
end

function out=normalize_coeff(in,n,nf,label)
% NORMALIZE_COEFF Execute the documented normalize_coeff operation.
%
% Syntax:
%   out=normalize_coeff(in,n,nf,label)
%
% Inputs:
%   in              : [documented value] Input required by the implemented function contract.
%   n               : [integer scalar or array] Discrete count or index required by the algorithm.
%   nf              : [integer scalar or array] Discrete count or index required by the algorithm.
%   label           : [documented value] Input required by the implemented function contract.
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
    if iscell(in), out=zeros(n,n,nf); for k=1:nf, out(:,:,k)=in{k}; end
    else, out=in; end
    if nf==1&&ismatrix(out), out=reshape(out,n,n,1); end
    if size(out,1)~=n||size(out,2)~=n||size(out,3)~=nf
        error('CRESTU:CoefficientShape','%s array has the wrong shape.',label);
    end
end
