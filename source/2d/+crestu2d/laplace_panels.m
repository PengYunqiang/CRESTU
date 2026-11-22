function [S,H,Sy,Sz,Hy,Hz]=laplace_panels(points,boundary)
%ROUND3_LAPLACE_PANELS Exact constant-panel log(r) direct-BIE operators.
% H integrates the SOURCE-normal derivative; normals point out of fluid.
% No boundary jump is inserted here. Smooth collocation adds pi*I separately.
    dy=points(:,1)-boundary.mid(:,1).'; dz=points(:,2)-boundary.mid(:,2).';
    ty=boundary.tangent(:,1).';tz=boundary.tangent(:,2).';
    ny=boundary.normal(:,1).';nz=boundary.normal(:,2).';
    x=dy.*ty+dz.*tz;a=dy.*ny+dz.*nz;
    half=boundary.length.'/2;up=half-x;um=-half-x;aa=abs(a);
    rp=up.^2+a.^2;rm=um.^2+a.^2;
    angle=atan2(up,aa)-atan2(um,aa);
    S=primitive(up,aa,rp)-primitive(um,aa,rm);
    H=-sign(a).*angle;
    if nargout>2
        assert(all(rp>0 & rm>0,'all'),'CRESTU2D:FiniteProbe','Derivative target lies on a panel endpoint.');
        sx=.5*log(rm./rp);sa=-H;
        hx=a.*(1./rp-1./rm);ha=up./rp-um./rm;
        Sy=sx.*ty+sa.*ny;Sz=sx.*tz+sa.*nz;
        Hy=hx.*ty+ha.*ny;Hz=hx.*tz+ha.*nz;
    end
end
function value=primitive(u,a,r2)
    value=.5*u.*log(r2)-u+a.*atan2(u,a);
    value(r2==0)=0;
end
