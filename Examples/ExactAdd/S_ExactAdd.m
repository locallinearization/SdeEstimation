function [m,P] = S_ExactAdd(t,t0,m0,P0,theta)
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

Cm = m0.*exp(0.5.*theta(1).*sqr(t0));
 m = Cm.*exp(-0.5.*theta(1).*sqr(t));

 q = 2.*theta(3)+1; 
Cp = (P0 - 0.5.*sqr(theta(4))./theta(1)) .*exp(theta(1).*sqr(t0)) - (sqr(theta(2)).*t0.^q)./q;
 P = exp(-theta(1).*sqr(t)) .* ( Cp + (sqr(theta(2)).*t.^q)./q ) + 0.5.*sqr(theta(4))./theta(1);
 

%------------------

function y=sqr(x)
y=x.*x;