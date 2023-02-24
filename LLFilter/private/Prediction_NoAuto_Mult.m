function [ms, Ps, ms1, Ps1] = Prediction_NoAuto_Mult(tt,tt0,m0,P0,theta,FileName,rtol,atol)
%
% This function computes the first two conditional moments of the equation:
%   dx = f(x,t)*dt + sum((B*x+b(t))*dw),
% where, f(x)=A*x+a(t), a(t)=a0*t+a1, b(t)=b0*t+b1, A and B are constant matrices, and
% a0, a1, b0 and b1 are vectors
% These moments satisface the equations:
%   dm/dt = A*m+a(t),  m(0)=m0
%   dP/dt = A*P + P*A' + a(t)*m' + m*a'(t)+ sum(B*P*B' + B*m'*b'(t) + b(t)*m*B + b(t)*b'(t)), P(0)=P0

%  References
%     "Simplified formulas for the mean and variance of linear stochastic differential equations",  
%     Jimenez J.C., Applied Mathematics Letters, Volume 49, November 2015, Pages 12-19. 

[f0,ge,A,B,a1,b1]=FileName(tt0,m0,theta);
h=tt-tt0;
b  = bxt(m0,ge,B);

d=size(m0,1);
d2=d*d;

I_d=eye(d);
bbT=zeros(d,d);
b1bT=zeros(d,d);
bb1T=zeros(d,d);
b1b1T=zeros(d,d);
BK{1}=zeros(d2,d2);
BK{2}=zeros(d2,d);
BK{3}=zeros(d2,d);
BK{4}=zeros(d2,d);
BK{5}=zeros(d2,d);
m=size(B,2);
for i=1:m 
   bbT = bbT + b{i}*b{i}';
   BK{1}= BK{1} + kron(B{i},B{i});
   BK{2}= BK{2} + kron(b{i},B{i});
   BK{3}= BK{3} + kron(B{i},b{i});
   
   BK{4}= BK{4} + kron(b1{i},B{i});
   BK{5}= BK{5} + kron(B{i},b1{i});
   b1bT = b1bT + b1{i}*b{i}';
   bb1T = bb1T + b{i}*b1{i}';
   b1b1T= b1b1T + b1{i}*b1{i}';
end
a0=f0-A*m0;
v0=zeros(d,1);
LT=[I_d v0 v0];
R=[v0; 0; 1];

AE = kron(I_d,A) + kron(A,I_d) + BK{1};

beta1 = vec(bbT);
beta2 = vec(b1bT + bb1T);
beta3 = vec(b1b1T);
beta4 = kron(I_d,a0) + kron(a0,I_d) + BK{2} + BK{3};
beta5 = kron(I_d,a1) + kron(a1,I_d) + BK{4} + BK{5};

beta1 = beta1 + beta4*m0;
beta2 = beta2 + beta5*m0;
beta4 = beta4*LT;
beta5 = beta5*LT;

clear a0 v0 LT bbT b1bT bb1T BK

LC=d+2;
C=[A a1 f0; zeros(2,LC)];
C(d+1,LC)=1;
M1 = [AE beta5 beta4 beta3 beta2 beta1];
M2 = [zeros(LC,d2) C eye(LC) zeros(LC,3)];
M3 = [zeros(LC,d2+LC) C zeros(LC,3)];
M4 = [zeros(1,d2+2*LC) 0 2 0];
M5 = [zeros(1,d2+2*LC) 0 0 1];
EM = h*[M1; M2; M3; M4; M5; zeros(1,d2+2*LC+3)]; 

clear AE beta5 beta4 beta3 beta2 beta1 C M1 M2 M3 M4 M5

% select p-p of Pade
tolr=[1.0e-6, 1.0e-9 1.0e-12 1.0e-15];
Table=[2 3 4 4; 3 4 5 6];   
nhM = norm(EM,'inf');
col = find(tolr>=rtol,1,'last');
if isempty(col)
   pd = 2;
else
   fil = (nhM>=1)+1;
   pd = Table(fil,col);
end
% scaling calculation
[~,e] = log2(nhM);
s = max(0,e+1);
% exponential calculation
EM= expm64v4(EM,pd,s);

Ps = EM(1:d2,1:d2)*vec(P0) + EM(1:d2,d2+LC+1:d2+2*LC)*R + EM(1:d2,end);
ms = m0 + EM(d2+1:d2+d,d2+d+2);
Ps = reshape(Ps,d,d);

if nargout>2
  EM  = EM*EM;
  Ps1 = EM(1:d2,1:d2)*vec(P0) + EM(1:d2,d2+d+3:d2+d+2+LC)*R + EM(1:d2,end);
  ms1 = m0 + EM(d2+1:d2+d,d2+d+2);
  Ps1 = reshape(Ps1,d,d);
end  


%----------------------------------

function v = vec(x)

[n,m]=size(x);
v=reshape(x,n*m,1);

