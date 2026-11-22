function [S,H,Sy,Sz,Hy,Hz]=hf_operators(points,b,depth,order)
%HF_OPERATORS D(top)/N(bottom) strip Green function, Delta G=2*pi*delta.
% G=log|tanh(a*(w-ws))|-log|tanh(a*(w-conj(ws)))|, a=pi/(4*h).
% Exact log panels and exact reflected panels; only a smooth remainder uses Gauss.
% All source normals point OUT OF FLUID. No boundary jump in these operators.
    gradients=nargout>2;
    if gradients,[S,H,Sy,Sz,Hy,Hz]=crestu2d.laplace_panels(points,b);
    else,[S,H]=crestu2d.laplace_panels(points,b);end
    mirror=b;mirror.mid(:,2)=-b.mid(:,2);mirror.tangent(:,2)=-b.tangent(:,2);mirror.normal(:,2)=-b.normal(:,2);
    if gradients
        [s,h,sy,sz,hy,hz]=crestu2d.laplace_panels(points,mirror);
        Sy=Sy-sy;Sz=Sz-sz;Hy=Hy-hy;Hz=Hz-hz;
    else,[s,h]=crestu2d.laplace_panels(points,mirror);end
    S=S-s;H=H-h;
    if isinf(depth),return;end
    alpha=pi/(4*depth);[t,wt]=crestu2d.gauss(order);
    w=complex(points(:,1),points(:,2));ns=complex(b.normal(:,1),b.normal(:,2)).';nm=conj(ns);
    for j=1:order
        src=b.mid+(t(j)/2*b.length).*b.tangent;ws=complex(src(:,1),src(:,2)).';
        u=w-ws;v=w-conj(ws);weight=(wt(j)/2*b.length).';
        [fu,du,ddu]=remainder(u,alpha);[fv,dv,ddv]=remainder(v,alpha);
        S=S+real(fu-fv).*weight;H=H+real(-du.*ns+dv.*nm).*weight;
        if gradients
            Sy=Sy+real(du-dv).*weight;Sz=Sz-imag(du-dv).*weight;
            Hy=Hy+real(-ddu.*ns+ddv.*nm).*weight;Hz=Hz+imag(ddu.*ns-ddv.*nm).*weight;
        end
    end
end
function [f,d,dd]=remainder(u,a)
    v=a*u;small=abs(v)<1e-3;f=zeros(size(u));d=f;dd=f;
    s=v(small);f(small)=-s.^2/3+7*s.^4/90-62*s.^6/2835;
    d(small)=a*(-2*s/3+14*s.^3/45-124*s.^5/945);
    dd(small)=a^2*(-2/3+14*s.^2/15-124*s.^4/189);
    s=v(~small);f(~small)=log(tanh(s)./s);
    d(~small)=a*(2./sinh(2*s)-1./s);
    dd(~small)=a^2*(-4*cosh(2*s)./sinh(2*s).^2+1./s.^2);
end
