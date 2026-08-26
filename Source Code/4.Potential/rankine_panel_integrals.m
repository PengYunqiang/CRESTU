function [G, dGdn] = rankine_panel_integrals(P_in, X, normal)
% RANKINE_PANEL_INTEGRALS Evaluate the analytic Rankine source and source-normal panel integrals.
%
% Syntax:
%   [G, dGdn] = rankine_panel_integrals(P_in, X, normal)
%
% Description:
%   The routine implements a component of the linear Rankine boundary-element formulation for incompressible, irrotational gravity-wave flow. Geometry, reflection parity, free-surface impedance, and complex phase follow the project convention exp(i*omega*t).
%
% Inputs:
%   P_in               - [4 x 3] or [3 x 4] Ordered panel vertices in global coordinates, [m].
%   X                  - [1 x 3] Field or collocation point in global coordinates, [m].
%   normal             - [1 x 3] Unit source-panel normal, dimensionless.
%
% Outputs:
%   G                  - [scalar] Surface integral of the Rankine kernel 1/r, [m].
%   dGdn               - [scalar] Source-normal derivative integral of the Rankine kernel, dimensionless.
%
% Governing Equations / Theory:
%   Green third identity, the Rankine kernel 1/r, linearized free-surface theory, and reflection symmetry as applicable.
%
% References:
%   - Newman, J. N. (1977), Marine Hydrodynamics; Hess and Smith (1964); project boundary-condition specification.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% --- 1. Validate Inputs and Initialize the Algorithm ---

    P = squeeze(P_in);
    if size(P, 1) == 3 && size(P, 2) == 4
        P = P'; % Normalize the vertex array to a 4-by-3 matrix.
    end

    p_center = mean(P, 1);

    e3 = normal / norm(normal);
    e1 = P(2, :) - P(1, :);
    e1 = e1 - dot(e1, e3) * e3;
    if norm(e1) < 1e-12
        e1 = P(3, :) - P(2, :);
        e1 = e1 - dot(e1, e3) * e3;
    end
    e1 = e1 / norm(e1);
    e2 = cross(e3, e1);

    T = [e1; e2; e3];

    P_loc = (P - p_center) * T';
    X_loc = (X - p_center) * T';

    x = X_loc(1);
    y = X_loc(2);
    z = X_loc(3);

    xi = [P_loc(:, 1); P_loc(1, 1)];
    eta = [P_loc(:, 2); P_loc(1, 2)];

    G = 0.0;
    dGdn = 0.0;
    eps_tol = 1e-12;

    for k = 1:4
        x1 = xi(k);     y1 = eta(k);
        x2 = xi(k + 1);   y2 = eta(k + 1);

        dx = x2 - x1;
        dy = y2 - y1;
        d12 = sqrt(dx^2 + dy^2);
        if d12 < eps_tol, continue; end

        r1 = sqrt((x - x1)^2 + (y - y1)^2 + z^2);
        r2 = sqrt((x - x2)^2 + (y - y2)^2 + z^2);

        cos_th = dx / d12;
        sin_th = dy / d12;

        t0 = (x1 - x) * sin_th - (y1 - y) * cos_th;

        denom_log = r1 + r2 - d12;
        if denom_log > eps_tol
            log_term = log((r1 + r2 + d12) / denom_log);
            G = G + t0 * log_term;
        end

    end

    % Robust oriented solid angle.  For source-normal differentiation,
    % d/dn_y(1/|X-y|) is minus the oriented angle seen from X.
    if abs(z) > eps_tol
        a = P(1, :) - X; b = P(2, :) - X; c = P(3, :) - X; d = P(4, :) - X;
        omega = triangle_solid_angle(a, b, c) + triangle_solid_angle(a, c, d);
        dGdn = -omega;
    end

    G = G - z * dGdn;

end

function omega = triangle_solid_angle(a, b, c)
% TRIANGLE_SOLID_ANGLE Evaluate the signed solid angle subtended by a triangular facet.
%
% Syntax:
%   omega = triangle_solid_angle(a, b, c)
%
% Description:
%   The routine implements a component of the linear Rankine boundary-element formulation for incompressible, irrotational gravity-wave flow. Geometry, reflection parity, free-surface impedance, and complex phase follow the project convention exp(i*omega*t).
%
% Inputs:
%   a                  - [1 x 3] First triangle vertex relative to the field point, [m].
%   b                  - [1 x 3] Second triangle vertex relative to the field point, [m].
%   c                  - [1 x 3] Third triangle vertex relative to the field point, [m].
%
% Outputs:
%   omega              - [scalar] Signed solid angle, [sr].
%
% Governing Equations / Theory:
%   Green third identity, the Rankine kernel 1/r, linearized free-surface theory, and reflection symmetry as applicable.
%
% References:
%   - Newman, J. N. (1977), Marine Hydrodynamics; Hess and Smith (1964); project boundary-condition specification.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% --- 1. Validate Inputs and Initialize the Algorithm ---

    numerator = dot(a, cross(b, c));
    denominator = norm(a) * norm(b) * norm(c) + dot(a, b) * norm(c) + ...
        dot(b, c) * norm(a) + dot(c, a) * norm(b);
    omega = 2 * atan2(numerator, denominator);
end
