function [phi,gy,gz]=hf_field(r,points)
%HF_FIELD Strict fluid points or exactly on the homogeneous flat boundaries.
    [S,H,Sy,Sz,Hy,Hz]=crestu2d.hf_operators(points,r.boundary,r.options.depth,r.options.quadrature_order);
    phi=(H*r.potential-S*r.normal_derivative)/(2*pi);
    gy=(Hy*r.potential-Sy*r.normal_derivative)/(2*pi);gz=(Hz*r.potential-Sz*r.normal_derivative)/(2*pi);
end
