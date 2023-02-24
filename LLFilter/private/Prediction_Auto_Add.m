function [ms, Ps, ms1, Ps1] = Prediction_Auto_Add(tt,tt0,m0,P0,theta,FileName,rtol,atol)
%
% This function computes the first two conditional moments of the equation:
%   dx = f(x)*dt + sum(b*dw),  
% where, f(x)=A*x+a, A is a constant matrix, and a and b are vectors
% These moments satisface the equations:
%   dm/dt = A*m+a,  m(0)=m0
%   dP/dt = A*P + P*A' + a*m' + m*a' + sum(b*b'), P(0)=P0

%  Reference
%     "Simplified formulas for the mean and variance of linear stochastic differential equations",  
%     Jimenez J.C., Applied Mathematics Letters, Volume 49, November 2015, Pages 12-19. 

[f0,b,A]  = FileName(tt0,m0,theta);
h=tt-tt0;

d=size(m0,1);

a0=f0-A*m0;
m=size(b,2);
bbT=zeros(d);
for i=1:m 
   bbT = bbT + b{i}*b{i}';
end
C=h*[A a0 a0*m0'+0.5.*bbT f0; zeros(1,d+1) f0' 0; zeros(d,d+1) -A' zeros(d,1); zeros(1,2*d+2)];

% select p-p of Pade
tolr=[1.0e-6, 1.0e-9 1.0e-12 1.0e-15];
Table=[2 3 4 4; 3 4 5 6];   
nhC = norm(C,'inf');
col = find(tolr>=rtol,1,'last');
if isempty(col)
   pd = 2;
else
   fil = (nhC>=1)+1;
   pd = Table(fil,col);
end
% scaling calculation
[~,e] = log2(nhC);
s = max(0,e+1);
% exponential calculation
M= expm64v4(C,pd,s);

F1=M(1:d,1:d);
g1=M(1:d,2*d+2);
H1=M(1:d,d+2:end-1);

Ps = H1*F1';
Ps = F1*P0*F1' + Ps + Ps';
ms = m0 + g1;

if nargout>2
  MM=M*M;
  F1=MM(1:d,1:d);
  g1=MM(1:d,2*d+2);
  H1=MM(1:d,d+2:end-1);

  Ps1 = H1*F1';
  Ps1 = F1*P0*F1' + Ps1 + Ps1';
  ms1 = m0 + g1;
end  

