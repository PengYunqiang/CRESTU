function [phi,gy,gz]=finite_field(result,points)
%FINITE_FIELD Direct representation at strictly interior fluid points.
% Input is one result.finite_state{frequency}, with all four saved columns.
% Availability does not certify local body/corner fields; see Gate B scope.
% No jump terms, no smoothing and no interpolation of the solved trace.
    if nargout>1
        [S,H,Sy,Sz,Hy,Hz]=crestu2d.finite_operators(points,result.boundary);
        gy=(Hy*result.potential-Sy*result.normal_derivative)/(2*pi);
        gz=(Hz*result.potential-Sz*result.normal_derivative)/(2*pi);
    else
        [S,H]=crestu2d.finite_operators(points,result.boundary);
    end
    phi=(H*result.potential-S*result.normal_derivative)/(2*pi);
end
