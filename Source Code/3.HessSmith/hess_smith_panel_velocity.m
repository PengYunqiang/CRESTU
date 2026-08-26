function [u, v, w] = hess_smith_panel_velocity(P, X, p_center, e1, e2, e3) %#ok<INUSD>
% HESS_SMITH_PANEL_VELOCITY Evaluate the velocity induced by a constant-strength planar source panel.
%
% Syntax:
%   [u, v, w] = hess_smith_panel_velocity(P, X, p_center, e1, e2, e3)
%
% Description:
%   The implementation uses the constant-strength Hess-Smith source-panel formulation in a panel-local orthonormal frame. Analytic edge logarithms and an oriented solid angle provide the tangential and normal influence components without changing the panel-source convention.
%
% Inputs:
%   P                  - [4 x 3] Ordered panel vertices in global Cartesian coordinates, [m].
%   X                  - [1 x 3] Field or collocation point in global coordinates, [m].
%   p_center           - [1 x 3] Panel centroid in global coordinates, [m].
%   e1                 - [1 x 3] First unit tangent vector, dimensionless.
%   e2                 - [1 x 3] Second unit tangent vector, dimensionless.
%   e3                 - [1 x 3] Unit panel-normal vector, dimensionless.
%
% Outputs:
%   u                  - [scalar] Global x component of the unit-source induced-velocity influence coefficient, dimensionless.
%   v                  - [scalar] Global y component of the unit-source induced-velocity influence coefficient, dimensionless.
%   w                  - [scalar] Global z component of the unit-source induced-velocity influence coefficient, dimensionless.
%
% Governing Equations / Theory:
%   Hess and Smith constant-source panel integration and oriented solid-angle geometry.
%
% References:
%   - Hess, J. L. and Smith, A. M. O. (1964), calculation of nonlifting potential flow about arbitrary three-dimensional bodies.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% --- 1. Validate Inputs and Initialize the Algorithm ---

    e3 = e3 / norm(e3);
    e1 = e1 - dot(e1, e3) * e3;
    if norm(e1) < 1e-14
        edge = P(2, :) - P(1, :);
        e1 = edge - dot(edge, e3) * e3;
    end
    e1 = e1 / norm(e1);
    e2 = cross(e3, e1);
    e2 = e2 / norm(e2);

    T = [e1; e2; e3];
    P_loc = (P - p_center) * T';
    X_loc = (X - p_center) * T';

    x = X_loc(1);
    y = X_loc(2);
    z = X_loc(3);

    % A constant flat panel is represented by the projection of the four
    % vertices onto the local panel plane.
    xi = [P_loc(:, 1); P_loc(1, 1)];
    eta = [P_loc(:, 2); P_loc(1, 2)];

    u_loc = 0.0;
    v_loc = 0.0;
    eps_tol = 1e-14;

    % In-plane components: standard Hess-Smith edge-logarithm formula.
    for k = 1:4
        x1 = xi(k);
        y1 = eta(k);
        x2 = xi(k + 1);
        y2 = eta(k + 1);

        dx = x2 - x1;
        dy = y2 - y1;
        edge_len = hypot(dx, dy);
        if edge_len < eps_tol
            continue;
        end

        r1 = sqrt((x - x1)^2 + (y - y1)^2 + z^2);
        r2 = sqrt((x - x2)^2 + (y - y2)^2 + z^2);

        denom = r1 + r2 - edge_len;
        numer = r1 + r2 + edge_len;
        if denom > eps_tol && numer > denom
            log_term = log(numer / denom);
            u_loc = u_loc + (dy / edge_len) * log_term;
            v_loc = v_loc - (dx / edge_len) * log_term;
        end
    end

    % Normal component: signed solid angle of two consistently oriented
    % triangles. This form has the correct far-field decay and avoids
    % inverse-tangent branch/complement errors.
    q1 = [P_loc(1, 1) - x, P_loc(1, 2) - y, -z];
    q2 = [P_loc(2, 1) - x, P_loc(2, 2) - y, -z];
    q3 = [P_loc(3, 1) - x, P_loc(3, 2) - y, -z];
    q4 = [P_loc(4, 1) - x, P_loc(4, 2) - y, -z];

    omega = triangle_solid_angle(q1, q2, q3) + ...
            triangle_solid_angle(q1, q3, q4);
    w_loc = -omega;

    vel_glob = ([u_loc, v_loc, w_loc] / (4.0 * pi)) * T;
    u = vel_glob(1);
    v = vel_glob(2);
    w = vel_glob(3);
end

function omega = triangle_solid_angle(a, b, c)
% TRIANGLE_SOLID_ANGLE Evaluate the signed solid angle subtended by a triangular facet.
%
% Syntax:
%   omega = triangle_solid_angle(a, b, c)
%
% Description:
%   The implementation uses the constant-strength Hess-Smith source-panel formulation in a panel-local orthonormal frame. Analytic edge logarithms and an oriented solid angle provide the tangential and normal influence components without changing the panel-source convention.
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
%   Hess and Smith constant-source panel integration and oriented solid-angle geometry.
%
% References:
%   - Hess, J. L. and Smith, A. M. O. (1964), calculation of nonlifting potential flow about arbitrary three-dimensional bodies.
%
% Lead Authors: Yunqiang Peng, Zhentao Jiang (SJTU)

%% --- 1. Validate Inputs and Initialize the Algorithm ---

    la = norm(a);
    lb = norm(b);
    lc = norm(c);

    numerator = dot(a, cross(b, c));
    denominator = la * lb * lc + ...
                  dot(a, b) * lc + ...
                  dot(b, c) * la + ...
                  dot(c, a) * lb;

    omega = 2.0 * atan2(numerator, denominator);
end
