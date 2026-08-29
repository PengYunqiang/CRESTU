function [cfg, report] = tune_sponge_layer(cfg, omegas)
% TUNE_SPONGE_LAYER Tune sponge layer for the CRESTU hydrodynamic workflow.
%
% Syntax:
%   [cfg, report] = tune_sponge_layer(cfg, omegas)
%
% Description:
%   Prepares boundary-panel geometry for the Rankine solver.
%   Global coordinates and panel-normal signs are preserved.
%
% Inputs:
%   cfg                - [struct] Validated CRESTU configuration containing SI-valued physical and numerical parameters.
%   omegas             - [1 x Nf] Angular-frequency samples, [rad/s].
%
% Outputs:
%   cfg                - [struct] Parsed and validated CRESTU configuration in SI units.
%   report             - [struct] Mesh-quality or tuning diagnostics with documented units.
%
% Governing Equations / Theory:
%   Planar panel geometry, polygon moments, reflection transformations, and mesh-topology relations as applicable.
%
% References:
%   - Hess and Smith (1964); CRESTU BMF mesh-format and geometry conventions.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    omegas = reshape(omegas, 1, []);
    wavelength = zeros(size(omegas));
    for k = 1:numel(omegas)
        wave_number = local_dispersion(omegas(k), cfg.grav, cfg.water_depth);
        wavelength(k) = 2 * pi / wave_number;
    end
    original_outer = cfg.fs.r_outer;
    required_width = 1.5 * max(wavelength);
    cfg.fs.r_outer = max(cfg.fs.r_outer, cfg.fs.r_inner + required_width);
    width = cfg.fs.r_outer - cfg.fs.r_inner;
% Quadratic profile: integral(mu dr)=mu0*width/3.  The following
% coefficient targets exp(-4) attenuation over the sponge.
    dynamic_mu0 = zeros(size(omegas));
    for k = 1:numel(omegas)
        wave_number = 2 * pi / wavelength(k);
        dynamic_mu0(k) = min(2.5, max(cfg.fs.mu0, 12 / (wave_number * width)));
    end
    report = struct('wavelengths', wavelength,'required_width', required_width, ...
'original_outer_radius', original_outer,'outer_radius', cfg.fs.r_outer, ...
'width', width,'dynamic_mu0', dynamic_mu0,'target_attenuation', exp(-4));
    if cfg.fs.r_outer > original_outer * (1 + 1e-12)
        fprintf('[WARN] Sponge radius expanded: Rout %.3g -> %.3g m (width %.3g m >= 1.5 lambda_max).\n', ...
            original_outer, cfg.fs.r_outer, width);
    end
end

function k = local_dispersion(omega, g, depth)
% LOCAL_DISPERSION Solve the local dispersion relation used for sponge-layer tuning.
%
% Syntax:
%   k = local_dispersion(omega, g, depth)
%
% Description:
%   Prepares boundary-panel geometry for the Rankine solver.
%   Global coordinates and panel-normal signs are preserved.
%
% Inputs:
%   omega              - [scalar] Angular frequency, [rad/s].
%   g                  - [scalar] Gravitational acceleration, [m/s^2].
%   depth              - [scalar] Water depth; a nonpositive value denotes infinite depth, [m].
%
% Outputs:
%   k                  - [scalar] Gravity-wave wavenumber, [1/m].
%
% Governing Equations / Theory:
%   Planar panel geometry, polygon moments, reflection transformations, and mesh-topology relations as applicable.
%
% References:
%   - Hess and Smith (1964); CRESTU BMF mesh-format and geometry conventions.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% Stage 1: Validate Inputs and Initialize the Algorithm

    if depth <= 0
        k = omega^2 / g;
        return;
    end
    k = max(omega^2 / g, omega / sqrt(g * depth));
    for iteration = 1:50
        kh = k * depth;
        th = tanh(kh);
        residual = g * k * th - omega^2;
        derivative = g * (th + kh / (cosh(kh)^2));
        update = residual / derivative;
        k = max(k - update, eps);
        if abs(update) < 1e-13 * max(1, k)
            break;
        end
    end
end
