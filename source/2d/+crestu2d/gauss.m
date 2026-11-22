function [x,w] = gauss(n)
%GAUSS Gauss-Legendre quadrature on [-1,1], Golub-Welsch construction.
    assert(isscalar(n) && n==fix(n) && n>=1,'CRESTU2D:Quadrature','Invalid quadrature order.');
    j=(1:n-1)'; b=j./sqrt(4*j.^2-1);
    [v,d]=eig(diag(b,1)+diag(b,-1));
    [x,ix]=sort(diag(d)); w=2*v(1,ix)'.^2;
end
