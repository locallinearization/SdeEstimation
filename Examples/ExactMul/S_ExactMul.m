function [m,P] = S_ExactMul(t,t0,m0,P0,theta)
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

A= theta(1);
B= theta(2);
t2=t*t;
t02=t0*t0;

m = m0.*exp(0.5.*A.*(t2-t02));
P = P0.*exp(0.5.*(2.*A+B*B).*(t2-t02));
 