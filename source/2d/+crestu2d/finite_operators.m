function [S,H,Sy,Sz,Hy,Hz]=finite_operators(points,boundary)
%FINITE_OPERATORS Exact Neumann-bottom image operators for the direct BEM.
% The source-normal mirror is (ny,-nz); its reflected tangent is (ty,-tz).
    if nargout>2
        [S,H,Sy,Sz,Hy,Hz]=crestu2d.laplace_panels(points,boundary);
    else
        [S,H]=crestu2d.laplace_panels(points,boundary);
    end
    assert(strcmp(boundary.bottom_representation,'image_neumann'),'CRESTU2D:FiniteBottom','Bottom image representation required.');
    if strcmp(boundary.bottom_representation,'image_neumann')
        mirror=boundary;
        mirror.mid(:,2)=-2*boundary.depth-boundary.mid(:,2);
        mirror.tangent(:,2)=-boundary.tangent(:,2);
        mirror.normal(:,2)=-boundary.normal(:,2);
        if nargout>2
            [s,h,sy,sz,hy,hz]=crestu2d.laplace_panels(points,mirror);
            Sy=Sy+sy;Sz=Sz+sz;Hy=Hy+hy;Hz=Hz+hz;
        else
            [s,h]=crestu2d.laplace_panels(points,mirror);
        end
        S=S+s;H=H+h;
    end
end
