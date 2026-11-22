function [e,rem] = wave_e(z)
%WAVE_E e=exp(-z)*Ei(z), principal Ei, Re(z)>0.
% Separate gamma+log(z) so the near-waterline logarithm is integrated exactly.
% Ei(z)=-E1(-z)+i*pi*sign(Im z) away from the real axis.
    gamma=0.5772156649015328606;
    assert(all(real(z)>0,'all'),'CRESTU2D:Kernel','Both Green-function points cannot lie on the waterline.');
    e=complex(zeros(size(z))); rem=e; small=abs(z)<0.5;
    zs=z(small); term=zs; series=zs;
    for j=2:24
        term=term.*zs/j; series=series+term/j;
    end
    logpart=gamma+log(zs);
    rem(small)=expm1(-zs).*logpart+exp(-zs).*series;
    e(small)=logpart+rem(small);
    zz=z(~small);
    ei=-expint(-zz)+1i*pi*sign(imag(zz));
    realaxis=imag(zz)==0; ei(realaxis)=real(ei(realaxis));
    e(~small)=exp(-zz).*ei;
    rem(~small)=e(~small)-gamma-log(zz);
    assert(all(isfinite(e),'all'),'CRESTU2D:Kernel','Nonfinite exponential integral; input outside supported numerical range.');
end
