function [ms, Ps, ms1, Ps1] = Prediction_ODE_K(tt,tt0,m0,P0,theta,FileName,rtol,atol)
%
% This function computes solution of the equation
%   dx = f(x)*dt , x(0)=x0
% where, f(x)=A*x+a, A is a constant matrix, and a is a vector

% Reference
%   "Locally Linearized Runge-Kutta method of Dormand and Prince 
%    for large systems of initial value problems"
%    F.S.Naranjo-Noda, J.C.Jimenez
%    Journal of Computational Physics 426 (2021) 109946 


persistent kdim

[f0,A]  = FileName(tt0,m0,theta);
h=tt-tt0;

if isempty(kdim)
  kdim = 2;
end
kdmin=2;
kdmax=length(m0);
[phi, kdim] = phi1_adapt(A,f0,h,m0,kdim,rtol,atol,kdmax,kdmin);

ms = m0 + phi(:,1);
Ps = ms*ms';

if nargout>2
  ms1 = m0 + phi(:,2);
  Ps1 = ms1*ms1';
end




