function comparison = compare_wamit_results(results,wamit_reference)
% COMPARE_WAMIT_RESULTS Execute the documented compare_wamit_results operation.
%
% Syntax:
%   comparison = compare_wamit_results(results,wamit_reference)
%
% Inputs:
%   results         : [documented value] Input required by the implemented function contract.
%   wamit_reference : [documented value] Input required by the implemented function contract.
%
% Outputs:
%   comparison      : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
%COMPARE_WAMIT_RESULTS Match frequencies and report A33/B33 errors/symmetry.
    nf=numel(results.omegas);
    template=struct('omega',0,'reference_omega',0,'A33',0,'A33_reference',0, ...
        'A33_relative_error',0,'B33',0,'B33_reference',0,'B33_relative_error',0, ...
        'A_symmetry_error',0,'B_symmetry_error',0);
    comparison=repmat(template,nf,1);
    for k=1:nf
        [~,q]=min(abs(wamit_reference.omegas-results.omegas(k)));
        A=results.added_mass(:,:,k); B=results.damping(:,:,k);
        Ar=wamit_reference.added_mass(:,:,q); Br=wamit_reference.damping(:,:,q);
        comparison(k)=struct('omega',results.omegas(k),'reference_omega',wamit_reference.omegas(q), ...
            'A33',A(3,3),'A33_reference',Ar(3,3),'A33_relative_error',relative_error(A(3,3),Ar(3,3)), ...
            'B33',B(3,3),'B33_reference',Br(3,3),'B33_relative_error',relative_error(B(3,3),Br(3,3)), ...
            'A_symmetry_error',norm(A-A.','fro')/max(norm(A,'fro'),eps), ...
            'B_symmetry_error',norm(B-B.','fro')/max(norm(B,'fro'),eps));
    end
end

function e=relative_error(value,reference)
% RELATIVE_ERROR Execute the documented relative_error operation.
%
% Syntax:
%   e=relative_error(value,reference)
%
% Inputs:
%   value           : [documented value] Input required by the implemented function contract.
%   reference       : [documented value] Input required by the implemented function contract.
%
% Outputs:
%   e               : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
    e=abs(value-reference)/max(abs(reference),eps);
end
