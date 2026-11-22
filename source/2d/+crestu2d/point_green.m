function [G,Gy,Gd] = point_green(k,y,d,ys,ds)
%POINT_GREEN Unscaled (+log r singularity) outgoing deep-water source.
    s=y-ys; v=d-ds; D=d+ds; r2=s.^2+v.^2; ri2=s.^2+D.^2;
    assert(k>0 && all(r2>0,'all'),'CRESTU2D:Kernel','Positive k and distinct field/source points required.');
    e=crestu2d.wave_e(k*(D+1i*s));
    wave=exp(-k*D); co=cos(k*s); si=sin(k*s);
    G=.5*log(r2./ri2)+2*real(e)+2i*pi*wave.*co;
    Gy=s./r2+s./ri2+2*k*imag(e)-2i*pi*k*wave.*si;
    Gd=v./r2+D./ri2-2*k*real(e)-2i*pi*k*wave.*co;
end
