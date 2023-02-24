function [ms, Ps, ms1, Ps1] = Prediction_ODE(tt,tt0,m0,P0,theta,FileName,rtol,atol)
%
% This function computes solution of the equation
%   dx = f(x)*dt , x(0)=x0
% where, f(x)=A*x+a, A is a constant matrix, and a is a vector

% Reference
%   Dynamic properties of the Local Linearization method for initial-value problems 
%   Jimenez J.C., Biscay R., Mora C.M., Rodriguez L.M.
%   Appl. Math. Comput., 126 (2002) 63-80. 

[f0,A]  = FileName(tt0,m0,theta);
h=tt-tt0;

d=size(m0,1);
C = [A f0; zeros(1,d+1)];

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
M = expm64v4(h*C,pd,s);

ms = m0 + M(1:d,end);
Ps = ms*ms';

if nargout>2
  M = M*M;
  ms1 = m0 + M(1:d,end);
  Ps1 = ms1*ms1';
end


