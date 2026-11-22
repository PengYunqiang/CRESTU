function S = lid_potential(target_y,endpoints_y,k,gauss_order)
%LID_POTENTIAL Exact-waterline source potential integrals in deep water.
% S(i,j)=int_{endpoints(j)}^{endpoints(j+1)} G(y_i,0;eta,0) d eta.
% The kernel is unscaled (+log r for a submerged source), exp(+i*omega*t).
% On the waterline G has 2*log|y-eta|. No artificial source/target depth.
% Constant/log/cos terms are analytic; the Ei remainder is split at ANY
% target inside a source panel, resolving its leading pi*k*abs(y-eta) cusp.
    assert(isnumeric(target_y) && isreal(target_y) && isvector(target_y) && ...
        ~isempty(target_y) && all(isfinite(target_y)), ...
        'CRESTU2D:LidKernel','Finite real target coordinates required.');
    assert(isnumeric(endpoints_y) && isreal(endpoints_y) && isvector(endpoints_y) && ...
        numel(endpoints_y)>=2 && all(isfinite(endpoints_y)) && all(diff(endpoints_y)>0), ...
        'CRESTU2D:LidKernel','Lid endpoints must be finite and strictly increasing.');
    assert(isnumeric(k) && isreal(k) && isscalar(k) && isfinite(k) && k>0 && ...
        isnumeric(gauss_order) && isreal(gauss_order) && isscalar(gauss_order) && ...
        isfinite(gauss_order) && gauss_order==fix(gauss_order) && gauss_order>=2, ...
        'CRESTU2D:LidKernel','Positive k and integer Gauss order >=2 required.');
    y=target_y(:);ends=endpoints_y(:).';a=ends(1:end-1);b=ends(2:end);
    len=b-a;mid=(a+b)/2;half=len/2;
    logIntegral=xlogabs(b-y)-xlogabs(a-y)-len;
    % Stable analytic log integral when the panel is remote from the target.
    ratio=(mid-y)./half;far=abs(ratio)>12;
    if any(far,'all')
        power=1./ratio(far).^2;series=zeros(size(power));term=power;
        for j=1:12
            series=series+term/(j*(2*j+1));term=term.*power;
        end
        widths=repmat(half,numel(y),1);delta=mid-y;
        logIntegral(far)=2*widths(far).*log(abs(delta(far)))-widths(far).*series;
    end
    [q,w]=crestu2d.gauss(gauss_order);
    split=min(max(y,a),b);
    h1=(split-a)/2;c1=(split+a)/2;
    h2=(b-split)/2;c2=(b+split)/2;
    remainderIntegral=zeros(numel(y),numel(a));
    % <<<CORE>>> R(Z)=exp(-Z)*Ei(Z)-gamma-log(Z), R(0)=0 by continuity.
    for j=1:numel(q)
        z1=1i*k*(y-(c1+h1*q(j)));
        z2=1i*k*(y-(c2+h2*q(j)));
        remainderIntegral=remainderIntegral+w(j)*(h1.*real(remainder(z1))+h2.*real(remainder(z2)));
    end
    phase=k*half;sinc0=ones(size(phase));nz=phase~=0;
    sinc0(nz)=sin(phase(nz))./phase(nz);
    cosineIntegral=(len.*sinc0).*cos(k*(y-mid));
    gamma=0.5772156649015328606;
    S=2*((gamma+log(k))*len+logIntegral+remainderIntegral)+2i*pi*cosineIntegral;
    % <<< /CORE >>>
    assert(all(isfinite(S),'all'),'CRESTU2D:LidKernel','Nonfinite exact-waterline integral.');
end

function value=xlogabs(x)
    value=zeros(size(x));nz=x~=0;value(nz)=x(nz).*log(abs(x(nz)));
end

function r=remainder(z)
% Principal Ei on the imaginary axis, as the limit from Re(z)>0.
% Ei(i*x)=-E1(-i*x)+i*pi*sign(x), for real nonzero x.
    gamma=0.5772156649015328606;r=complex(zeros(size(z)));
    small=abs(z)<.5 & z~=0;zs=z(small);
    term=zs;series=zs;
    for j=2:24
        term=term.*zs/j;series=series+term/j;
    end
    logpart=gamma+log(zs);
    r(small)=expm1(-zs).*logpart+exp(-zs).*series;
    large=~small & z~=0;zl=z(large);
    ei=-expint(-zl)+1i*pi*sign(imag(zl));
    r(large)=exp(-zl).*ei-gamma-log(zl);
    % Entries with z==0 retain the analytic remainder limit 0.
end
