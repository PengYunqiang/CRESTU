function reference = read_wamit_excitation(filename,rho,g,wave_amplitude,characteristic_length)
% READ_WAMIT_EXCITATION Execute the documented read_wamit_excitation operation.
%
% Syntax:
%   reference = read_wamit_excitation(filename,rho,g,wave_amplitude,characteristic_length)
%
% Inputs:
%   filename        : [char|string] Input or output file path.
%   rho             : [scalar] Water density, in kg/m^3.
%   g               : [scalar] Gravitational acceleration, in m/s^2.
%   wave_amplitude  : [scalar] Incident-wave amplitude, in m.
%   characteristic_length: [documented value] Input required by the implemented function contract.
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
%READ_WAMIT_EXCITATION Read WAMIT .2 exciting-force coefficients in SI.
% Translational normalization is rho*g*A*L^2; moments use one more L.
    if nargin<4||isempty(wave_amplitude), wave_amplitude=1; end
    if nargin<5||isempty(characteristic_length), characteristic_length=1; end
    raw=readmatrix(filename,'FileType','text');
    if size(raw,2)<7, error('CRESTU:WamitFormat','Expected seven columns in %s.',filename); end
    periods=unique(raw(:,1),'stable'); headings=unique(raw(:,2),'stable');
    nf=numel(periods); nh=numel(headings); ndof=max(raw(:,3));
    value=complex(zeros(ndof,nh,nf)); coefficient=complex(zeros(ndof,nh,nf));
    for r=1:size(raw,1)
        k=find(periods==raw(r,1),1); h=find(headings==raw(r,2),1); mode=raw(r,3);
        coefficient(mode,h,k)=complex(raw(r,6),raw(r,7));
        scale=rho*g*wave_amplitude*characteristic_length^(2+(mode>3));
        value(mode,h,k)=scale*coefficient(mode,h,k);
    end
    reference=struct('file',filename,'periods',periods(:).','omegas',2*pi./periods(:).', ...
        'headings',headings(:).','coefficient',coefficient,'force',value);
end
