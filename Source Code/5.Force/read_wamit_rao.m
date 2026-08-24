function reference = read_wamit_rao(filename)
% READ_WAMIT_RAO Execute the documented read_wamit_rao operation.
%
% Syntax:
%   reference = read_wamit_rao(filename)
%
% Inputs:
%   filename        : [char|string] Input or output file path.
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
%READ_WAMIT_RAO Read WAMIT .4 complex motion RAOs.
    raw=readmatrix(filename,'FileType','text');
    if size(raw,2)<7, error('CRESTU:WamitFormat','Expected seven columns in %s.',filename); end
    periods=unique(raw(:,1),'stable'); headings=unique(raw(:,2),'stable');
    nf=numel(periods); nh=numel(headings); ndof=max(raw(:,3)); value=complex(zeros(ndof,nh,nf));
    for r=1:size(raw,1)
        k=find(periods==raw(r,1),1); h=find(headings==raw(r,2),1); mode=raw(r,3);
        value(mode,h,k)=complex(raw(r,6),raw(r,7));
    end
    reference=struct('file',filename,'periods',periods(:).','omegas',2*pi./periods(:).', ...
        'headings',headings(:).','complex',value,'amplitude',abs(value),'phase_deg',angle(value)*180/pi);
end
