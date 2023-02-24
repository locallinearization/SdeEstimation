function error=error_step(y0,ynew,y1,d,Atol,Rtol)
l=max(abs(y0),abs(ynew));
sc=Atol+Rtol.*l;
v=(ynew-y1)./sc;
error=sqrt((1./d).*(v'*v));
