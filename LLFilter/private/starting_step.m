function hinic=starting_step(f0,fx0,ft0,y0,d,Atol,Rtol,pp)
% pp=1/(1+order)   with order= 1 or 2 

q=1./d;
sc=Atol+Rtol.*abs(y0);
v0=y0./sc;
d0=sqrt(q.*(v0'*v0));
v1=f0./sc;
d1=sqrt(q.*(v1'*v1));
r=ft0+fx0*f0;
v2=r./sc;
d2=sqrt(q.*(v2'*v2));
temp=10.*Atol;
if d0<temp || d1<temp
    h0=Atol;
else h0=1.0e-2.*(d0./d1);
end
m=max(d1,d2);
if m<=eps
    h1=max(Atol,h0.*Rtol);
else h1=(0.01./m).^pp;
end
hinic=min(100.*h0,h1);