function reference = read_wamit_first_order(file_one,rho)
% READ_WAMIT_FIRST_ORDER Execute the documented read_wamit_first_order operation.
%
% Syntax:
%   reference = read_wamit_first_order(file_one,rho)
%
% Inputs:
%   file_one        : [documented value] Input required by the implemented function contract.
%   rho             : [scalar] Water density, in kg/m^3.
%
% Outputs:
%   reference       : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
%READ_WAMIT_FIRST_ORDER Read WAMIT .1 added-mass/damping coefficient data.
% For IPERIO=1 period output, col1=T. WAMIT coefficients are converted from
% A/rho and B/(rho*omega) to SI kg and kg/s (rotational units analogous).
    raw=readmatrix(file_one,'FileType','text');
    if size(raw,2)<5, error('CRESTU:WamitFormat','Expected five columns in %s.',file_one); end
    periods=unique(raw(:,1),'stable'); nf=numel(periods); ndof=max(max(raw(:,2:3)));
    A=zeros(ndof,ndof,nf); B=zeros(ndof,ndof,nf); omegas=2*pi./periods(:).';
    for r=1:size(raw,1)
        k=find(periods==raw(r,1),1); i=raw(r,2); j=raw(r,3);
        A(i,j,k)=rho*raw(r,4); B(i,j,k)=rho*omegas(k)*raw(r,5);
    end
    reference=struct('file',file_one,'periods',periods(:).','omegas',omegas, ...
        'added_mass',A,'damping',B,'rho',rho);
end
