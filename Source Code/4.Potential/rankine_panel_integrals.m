function [G, dGdn] = rankine_panel_integrals(P_in, X, normal)
% RANKINE_PANEL_INTEGRALS Execute the documented rankine_panel_integrals operation.
%
% Syntax:
%   [G, dGdn] = rankine_panel_integrals(P_in, X, normal)
%
% Inputs:
%   P_in            : [documented value] Input required by the implemented function contract.
%   X               : [documented value] Input required by the implemented function contract.
%   normal          : [integer scalar or array] Discrete count or index required by the algorithm.
%
% Outputs:
%   G               : [scalar] Gravitational acceleration, in m/s^2.
%   dGdn            : [documented value] Function result; dimensions and units follow the implemented contract.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
%
%
%

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
    
    xi  = [P_loc(:, 1); P_loc(1, 1)];
    eta = [P_loc(:, 2); P_loc(1, 2)];
    
    G = 0.0;
    dGdn = 0.0;
    eps_tol = 1e-12;
    
    for k = 1:4
        x1 = xi(k);     y1 = eta(k);
        x2 = xi(k+1);   y2 = eta(k+1);
        
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
        a=P(1,:)-X; b=P(2,:)-X; c=P(3,:)-X; d=P(4,:)-X;
        omega=triangle_solid_angle(a,b,c)+triangle_solid_angle(a,c,d);
        dGdn=-omega;
    end
    
    G = G - z * dGdn;

end

function omega=triangle_solid_angle(a,b,c)
% TRIANGLE_SOLID_ANGLE Execute the documented triangle_solid_angle operation.
%
% Syntax:
%   omega=triangle_solid_angle(a,b,c)
%
% Inputs:
%   a               : [documented value] Input required by the implemented function contract.
%   b               : [documented value] Input required by the implemented function contract.
%   c               : [documented value] Input required by the implemented function contract.
%
% Outputs:
%   omega           : [scalar] Angular frequency, in rad/s.
%
% Mathematical Reference:
%   See the inline equations and the corresponding CRESTU module theory notes.
%
% ==========================================
% Function implementation
% ==========================================
    numerator=dot(a,cross(b,c));
    denominator=norm(a)*norm(b)*norm(c)+dot(a,b)*norm(c)+ ...
        dot(b,c)*norm(a)+dot(c,a)*norm(b);
    omega=2*atan2(numerator,denominator);
end
