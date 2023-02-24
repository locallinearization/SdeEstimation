function [ms, Ps, ms1, Ps1] = Prediction_Auto_Mult_K(tt,tt0,m0,P0,theta,FileName,rtol,atol)
%
% This function computes the first two conditional moments of the equation:
%   dx = f(x)*dt + sum((B*x+b)*dw),
% where, f(x)=A*x+a, A and B are constant matrices, and a and b are constant vectors
% These moments satisface the equations:
%   dm/dt = A*m+a,  m(0)=m0
%   dP/dt = A*P + P*A' + a*m' + m*a' + sum(B*P*B' + B*m'*b' + b*m*B + b*b'), P(0)=P0

%  References
%     - "Simplified formulas for the mean and variance of linear stochastic differential equations",  
%       Jimenez J.C., Appl. Math. Letters, Volume 49, November 2015, Pages 12-19. 
%     - "Computing high dimensional multiple integrals involving matrix exponentials",  
%       Naranjo-Noda F.S. and Jimenez J.C., J. Comput. Appl. Math., 421 (2023) 114844.

persistent kdim

[f0,ge,A,B] = FileName(tt0,m0,theta);
h=tt-tt0;
b  = bxt(m0,ge,B);
 
d=size(m0,1);
d2=d*d;

I_d=eye(d);
bbT=zeros(d,d);
BK{1}=zeros(d2,d2);
BK{2}=zeros(d2,d);
BK{3}=zeros(d2,d);
m=max([size(B,2) size(b,2)]);
for i=1:m 
   bbT = bbT + b{i}*b{i}';
   BK{1}= BK{1} + kron(B{i},B{i});
   BK{2}= BK{2} + kron(b{i},B{i});
   BK{3}= BK{3} + kron(B{i},b{i}); 
end
a0=f0-A*m0;
v0=zeros(d,1);
LT=[I_d v0];

AE = kron(I_d,A) + kron(A,I_d) + BK{1};

beta1 = vec(bbT);
beta4 = kron(I_d,a0) + kron(a0,I_d) + BK{2} + BK{3};

beta1 = beta1 + beta4*m0;
beta4 = beta4*LT;

C = [A f0; zeros(1,d+1)];
M = [AE beta1 beta4; zeros(1,d2+d+2); zeros(d+1,d2+1) C];
r = [vec(P0); 1; zeros(d,1); 1];

if isempty(kdim)
  kdim = 2;
end
kdmin=2;
kdmax=length(r);
[phi, kdim] = phi1_adapt(M,M*r,h,r,kdim,rtol,atol,kdmax,kdmin);
phi = [r r] + phi;

Ps = phi(1:d2,1);
ms = m0 + phi(d2+2:d2+1+d,1);
Ps = reshape(Ps,d,d);

if nargout>2
  Ps1 = phi(1:d2,2);
  ms1 = m0 + phi(d2+2:d2+1+d,2);
  Ps1 = reshape(Ps1,d,d);
end

%----------------------------------

function v = vec(x)

[n,m]=size(x);
v=reshape(x,n*m,1);
