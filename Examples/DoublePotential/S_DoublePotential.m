function [m,P] = S_Ozaki(t,t0,m0,P0,theta)
% 
% Definition of the exact conditional predictions of the state variables
%
%  Input
%       t : final time                     (1 x 1)
%       t0: initial time                   (1 x 1)
%       m0: first moment at t0             (n_v_e x 1)
%       P0: second moment at t0            (n_v_e x n_v_e)
%    theta: model parameters               (1 x n_theta)
%  
%  Output
%       m: first moment at t               (n_v_e x 1)
%       P: second moment at t              (n_v_e x n_v_e)

fname = @DoublePotential;
options  = odeset('RelTol',1.0e-12,'AbsTol',1.0e-12);  
P0 = P0 - m0*m0';                          % to get the variance of the state variable
[tt,y] = ode15s(fname,[t0 t],[m0 P0],options,theta);
yy=y(end,:);   
m=yy(1);
P=yy(2) + m*m';                             % to return the 2nd moment of the state variable


function f=DoublePotential(t,x,a)
% expressions for the mean and variance of the state variable (NOT for the 2nd moment)

f=zeros(2,1);

f(1) = - a(1)*x(1) - a(2)*(3*x(2)*x(1) + x(1)^3);
f(2) = - 2*a(1)*x(2) - 6*a(2)*(x(2)^2 + x(2)*x(1)^2) + a(3)^2; 
