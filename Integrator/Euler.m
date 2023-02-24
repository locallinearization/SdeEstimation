function [t,z,w]=Euler(t0,x0,n_m,deltat,FileNames,theta,seed)
%   Euler-Maruyama scheme for the numerical integrator of the stochastic differential
% equations with multiplicative noise: 
%	  dx = f(t,x,theta) dt + g(t,x,theta) dw
% on uniform time partitions
%
% Notation
%   n_v_e: number of state variables x 
%
% Imputs:
%   t0:     initial time;
%   x0:     initial value of x at t0  (n_v_e x 1) 
%   n_m:    number of time-instants to integrate
%   deltat: integration step-size
%   FileNames: cell variable with the name of all the functions that define SS Model
%   theta:  parameters of the differential equation
%   seed:   seed for initializing the random number generator
%
% Outputs:
%   t:   integrations times (1 x n_m)
%   z:   numerical solucion obtained by the LL scheme (n_v_e x n_m)

ModEq = FileNames{1};

[~,ge] = ModEq(t0,x0,theta);
n_w = length(ge);

if (nargin==7) & (~isempty(seed))
   randn('seed',seed);
end
eta = randn(n_m,n_w);    % system noise generation 

n_v_e = size(x0,1);
t = zeros(1,n_m);
z = zeros(n_v_e,n_m);

x = x0;
te = t0;

w=sqrt(deltat) .* cumsum([zeros(1,n_w);eta])';
dw=diff(w,1,2);

t(1) = te;
z(:,1) = x;
for tt = 1:n_m-1
   [fe,ge] = ModEq(te,x,theta);
     
   m=size(ge,2);
   ruido=zeros(n_v_e,1);
   for k=1:m 
      ruido=ruido + ge{k}.*dw(k,tt);
   end 
   te = tt .* deltat + t0;
   x = x + fe.*deltat + ruido;

   t(tt+1) = te;
   z(:,tt+1) = x;
   
end% for en tt

