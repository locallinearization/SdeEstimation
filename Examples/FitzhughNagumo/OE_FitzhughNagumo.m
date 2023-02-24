function [he,Jhe,See] = OE_FitzhughNagumo(t,x,theta)
%
% Definition of the dicrete observation equation 
%            [he,Jhe,See] = OE(t,x,theta)  for observations with additive noise
% [he,Jhe,See,pe,Jpe,Lee] = OE(t,x,theta)  for observations with additive and multiplicative noise
%
%  Input 
%        t: time instant                                        (1 x 1)         
%        x: state variable at the time t                        (n_v_e x 1)   
%    theta: model parameters                                    (1 x n_theta)
%
%  Output 
%      he = h(t,x,theta)                                        (n_obs x 1)
%     Jhe : jacobian matrix of he                               (n_obs x n_obs)
%     See : variance of the additive observation noise          (n_obs x n_obs)
%      pe = p(t,x,theta)                                         cell(1,n) of (n_obs x 1)
%     Jpe : jacobian matrix of pe                                cell(1,n) of (n_obs x n_obs)
%     Lee : variance of the multiplicative observation noises   (1 x n) 

 he = x;   
Jhe = [1 0; 0 1];   

See=[theta(4) 0; 0 theta(5)];
