function [cfg,report] = tune_sponge_layer(cfg,omegas)
% TUNE_SPONGE_LAYER Execute the documented tune_sponge_layer operation.
%
% Syntax:
%   [cfg,report] = tune_sponge_layer(cfg,omegas)
%
% Inputs:
%   cfg             : [struct] Validated CRESTU configuration with SI-valued physical fields.
%   omegas          : [1 x Nf] Angular frequencies, in rad/s.
%
% Outputs:
%   cfg             : [struct] Validated CRESTU configuration with SI-valued physical fields.
%   report          : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
%TUNE_SPONGE_LAYER Enforce >=1.5 wavelengths and prepare frequency damping.
    omegas=reshape(omegas,1,[]); wavelength=zeros(size(omegas));
    for k=1:numel(omegas)
        wave_number=local_dispersion(omegas(k),cfg.grav,cfg.water_depth);
        wavelength(k)=2*pi/wave_number;
    end
    original_outer=cfg.fs.r_outer;
    required_width=1.5*max(wavelength);
    cfg.fs.r_outer=max(cfg.fs.r_outer,cfg.fs.r_inner+required_width);
    width=cfg.fs.r_outer-cfg.fs.r_inner;
    % Quadratic profile: integral(mu dr)=mu0*width/3.  The following
    % coefficient targets exp(-4) attenuation over the sponge.
    dynamic_mu0=zeros(size(omegas));
    for k=1:numel(omegas)
        wave_number=2*pi/wavelength(k);
        dynamic_mu0(k)=min(2.5,max(cfg.fs.mu0,12/(wave_number*width)));
    end
    report=struct('wavelengths',wavelength,'required_width',required_width, ...
        'original_outer_radius',original_outer,'outer_radius',cfg.fs.r_outer, ...
        'width',width,'dynamic_mu0',dynamic_mu0,'target_attenuation',exp(-4));
    if cfg.fs.r_outer>original_outer*(1+1e-12)
        fprintf('>>> Sponge radius expanded: Rout %.3g -> %.3g m (width %.3g m >= 1.5 lambda_max).\n', ...
            original_outer,cfg.fs.r_outer,width);
    end
end

function k=local_dispersion(omega,g,depth)
% LOCAL_DISPERSION Execute the documented local_dispersion operation.
%
% Syntax:
%   k=local_dispersion(omega,g,depth)
%
% Inputs:
%   omega           : [scalar] Angular frequency, in rad/s.
%   g               : [scalar] Gravitational acceleration, in m/s^2.
%   depth           : [scalar] Water depth, in m; zero denotes infinite depth.
%
% Outputs:
%   k               : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
    if depth<=0, k=omega^2/g; return; end
    k=max(omega^2/g,omega/sqrt(g*depth));
    for iteration=1:50
        kh=k*depth; th=tanh(kh); residual=g*k*th-omega^2;
        derivative=g*(th+kh/(cosh(kh)^2)); update=residual/derivative;
        k=max(k-update,eps);
        if abs(update)<1e-13*max(1,k), break; end
    end
end
