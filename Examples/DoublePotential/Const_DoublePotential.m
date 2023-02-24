% Script for model definition: nonlinear SDE with bimodal stationary distribution
%   H. Singer, Parameter estimation of nonlinear stochastic differential equations: 
%   Simulated maximum likelihood versus extended Kalman filter and Ito-Taylor expansion, 
%   J. Comput. Graph. Stats., 11 (2002) 972-995.

deltat_min=0.001;   % sampling period for integration
deltat=1;          % sampling period for filtering
t0=0;              % initial time 
tTotal=500+t0;     % sampling time  

% definition of state space parameters
Theta(1)=-1; 
Theta(2)=0.1; 
Theta(3)=2; 
Theta(4)=0;        % observation noise variance   

% index of the observed state variables
ObsX=[true];

% initial state variable for integration;
X0=0;
% initial filter estimates
Xf0=X0;      
% initial filter variance estimates
Pf0=0;         

% definition of the model equation 
ModEq = 'SE_DoublePotential';
% definition of the observation equation
ObsEq = 'OE_DoublePotential';
% definition of the exact predictions of the model
ExactSol = 'S_DoublePotential';

FileNames=MakeFileNames(ModEq,ObsEq,ExactSol);

% absolute and relative tolerances for filtering
AbsTol=[1.0e-6,1.0e-6,0.001];
% absolute and relative tolerances for filtering
RelTol=[1.0e-3,1.0e-3,0.001];




