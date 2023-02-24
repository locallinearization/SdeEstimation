function   [f,g,Jf,Jg] = SE_CIR(t,x,theta)
%
% Definition of the state equation
%               [f,Jf] = SE(t,x,theta)   for autonomous ODEs
%             [f,g,Jf] = SE(t,x,theta)   for autonomous SDEs with additive noise 
%          [f,g,Jf,Jg] = SE(t,x,theta)   for autonomous SDEs with multplicative noise
%       [f,g,Jf,ft,gt] = SE(t,x,theta)   for non-autonomous SDEs with additive noise
%    [f,g,Jf,Jg,ft,gt] = SE(t,x,theta)   for non-autonomous SDEs with multplicative noise
% [[f;1],[Jf ft; 0 0]] = SE(t,x,theta)   for non-autonomous ODEs
% 
%  Input 
%        t: time instant                             (1 x 1)         
%        x: state variable at the time t             (n_v_e x 1)   
%    theta: model parameters                         (1 x n_theta)
%
%  Output 
%       f: drift at t,x                              (n_v_e x 1)
%      Jf: jacobian of f                             (n_v_e x n_v_e)
%       g: difussion at t,x                          cell(1,m) of (n_v_e x 1)
%      Jg: jacobian of g                             cell(1,m) of (n_v_e x n_v_e)
%      ft: derivative with respect the time at t,x   (n_v_e x 1)
%      gt: derivative with respect the time at t,x   cell(1,m) of (n_v_e x 1)


f(1)=theta(1)+theta(2).*x;
Jf=theta(2);

g{1}=theta(3).*sqrt(x);
Jg{1}= 0.5.*theta(3)./sqrt(x);
