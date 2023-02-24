function [m,P] = S_CIR(t,t0,m0,P0,theta)
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

m = m0.*exp(theta(2).*(t-t0)) + (exp(theta(2).*(t-t0))-1).*(theta(1)./theta(2));

cte=(2.*theta(1)+theta(3).^2)./(2.*theta(2).^2);
P = cte.*((theta(1)+2.*m0.*theta(2)+P0./cte).*exp(2.*theta(2)*(t-t0))-2.*(theta(1)+m0.*theta(2)).*exp(theta(2).*(t-t0))+theta(1));