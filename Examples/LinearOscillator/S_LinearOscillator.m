function [m,P] = S_LinearOscillator(t,t0,m0,P0,theta)
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

a=theta(1);
w=theta(2);
s1=theta(3);
s2=theta(4);

M=[    0  1  1 0               0  0   0         0;...
      -w  0  0 1         a*m0(1)  a   0         0;...
      -w  0  0 1         a*m0(1)  a   0         0;...
   s1*s1 -w -w 0 2*a*m0(2)+s2*s2  0 2*a         0;...
       0  0  0 0               0  0   0         0;...
       0  0  0 0               0  0   1     m0(2);...
       0  0  0 0               0 -w   0 a-w*m0(1);...
       0  0  0 0               0  0   0         0];

u=[P0(1,1); P0(1,2); P0(1,2); P0(2,2); 1; 0; 0; 1];   

v=expm(M*(t-t0))*u;

d=length(m0);
d2=d*d;
L2=[zeros(d,d2+1) eye(d) zeros(d,1)];
L1=[eye(d2) zeros(d2,d+2)];

m=m0+L2*v;
P=reshape(L1*v,d,d);
 
       
   