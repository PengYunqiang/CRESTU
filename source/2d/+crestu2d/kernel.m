function [S,Sy,Sd] = kernel(target,geo,k,order,domain,own,tlocal)
%KERNEL Exact straight-panel logarithms and PV inverse; Gauss wave remainder.
% Entries integrate G and its field derivatives over physical source arclength.
% own identifies the panel containing each target; zero means off-body.
    if nargin<6, own=zeros(numel(target),1); tlocal=zeros(numel(target),1); end
    target=target(:); b=geo.b.'; c=geo.c.';
    [ld,id]=line_integrals(target-c,-b);
    ii=find(own>0);
    if ~isempty(ii)
        ix=sub2ind(size(ld),ii,own(ii));
        tt=tlocal(ii); bb=geo.b(own(ii));
        ld(ix)=2*log(abs(bb))+(1-tt).*log(1-tt)+(1+tt).*log(1+tt)-2;
        id(ix)=log((1+tt)./(1-tt))./bb;
    end
    if strcmp(domain,'unbounded')
        S=ld.*geo.jac.'; Sy=-imag(id).*geo.jac.'; Sd=real(id).*geo.jac.'; return;
    end
    [li,iv]=line_integrals(target+conj(c),conj(b));
    [q,weights]=crestu2d.gauss(order);
    rr=zeros(size(ld)); angleint=rr; wav=complex(rr);
    % <<<CORE>>> Green remainder; no artificial free-surface/outer boundary.
    for j=1:order
        z=target+conj(c+b*q(j));
        [~,rem]=crestu2d.wave_e(k*z);
        rr=rr+weights(j)*real(rem);
        angleint=angleint+weights(j)*(angle(z)+imag(rem));
        wav=wav+weights(j)*exp(-k*z);
    end
    constant=0.5772156649015328606+log(k);
    S=ld+li+4*constant+2*rr+2i*pi*real(wav);
    Sy=-imag(id)-imag(iv)+2*k*angleint+2i*pi*k*imag(wav);
    Sd=real(id)+real(iv)-2*k*(2*constant+li+rr)-2i*pi*k*real(wav);
    S=S.*geo.jac.'; Sy=Sy.*geo.jac.'; Sd=Sd.*geo.jac.';
    % <<< /CORE >>>
end
function [l,iv]=line_integrals(a,b)
% Integrals of log|a+b*t| and 1/(a+b*t), -1<t<1.
    q=a./b; qp=q+1; qm=q-1;
    l=2*log(abs(b))+real(xlog(qp)-xlog(qm)-2);
    iv=(log(qp)-log(qm))./b;
    % In the remote-panel limit avoid subtraction of large primitive values.
    far=abs(q)>12;
    if any(far,'all')
        qq=q(far); power=1./qq.^2; series=zeros(size(qq)); invseries=1./qq;
        for j=1:12
            series=series+power/(j*(2*j+1));
            invseries=invseries+power./qq/(2*j+1);
            power=power./qq.^2;
        end
        bb=repmat(b,size(a,1),1);
        l(far)=2*log(abs(a(far)))-real(series);
        iv(far)=2*invseries./bb(far);
    end
end
function y=xlog(x), y=x.*log(x); y(x==0)=0; end
