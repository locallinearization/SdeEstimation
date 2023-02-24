% Script for model definition: Cox-Ingersoll-Ross (CIR) model of interest rate. 
%   Cox, J. C., Ingersoll, J. E., and Ross, S. A. (1985) A theory of the term structure of interest rates,
%   Econometrica 53, 285–408.
% Model setting from: 
%   Chan, K. C., at el. An empirical comparison of alternative models of the short-term 
%   interest rate, J. Finance 47 (1992) 1209–1227.

deltat_min=0.00083;    % sampling period for integration
deltat=0.083;          % sampling period for filtering
t0=0;                  % initial time  
tTotal=t0+25.398;      % sampling time  

% definition of state space parameters
Theta(1)= 0.0189;
Theta(2)=-0.2339;
Theta(3)=sqrt(0.0073);
Theta(4)=0;              % observation noise variance

% index of the observed state variables
ObsX=[true];

% initial state variable for integration
X0=0.1;            
% initial filter estimates
Xf0=X0;                
% initial filter variance estimates
Pf0=0.0001;                   

% definition of the model equation 
ModEq = 'SE_CIR';
% definition of the observation equation
ObsEq = 'OE_CIR';
% definition of the exact predictions of the model
ExactSol = 'S_CIR';

FileNames=MakeFileNames(ModEq,ObsEq,ExactSol);

% absolute and relative tolerances for filtering
AbsTol=[1e-6,1e-6,0.001];
RelTol=[1e-3,1e-3,0.001];




