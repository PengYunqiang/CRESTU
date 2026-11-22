function nodes=hf_rectangle(b,d,n)
%HF_RECTANGLE Cosine-graded straight edges; exact corners, no geometry rounding.
    t=(1-cos(linspace(0,pi,n+1)'))/2;
    nodes=[-b/2*ones(n+1,1),-d*t; -b/2+b*t(2:end),-d*ones(n,1); b/2*ones(n,1),-d+d*t(2:end)];
    nodes(end,2)=0;
end
