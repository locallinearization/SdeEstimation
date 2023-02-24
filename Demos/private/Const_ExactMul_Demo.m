% Script for model definition: Non-autonomous linear SDE with multiplicative noise 
%     "Bias reduction in the estimation of diffusion processes from discrete observations",
%      Jimenez J.C, IMA Journal of Mathematical Control and Information, 
%      Volumen 37, December 2020, Pages 1468–1505. 

deltat_min=0.0001;  % sampling period for integration
deltat=1;           % sampling period for filtering
t0=0.5;             % initial time  
tTotal=10+t0;       % sampling time  

% definition of state space parameters
Theta(1)=-0.1;
Theta(2)=0.1;
Theta(3)=0.0001;    % observation noise variance

% index of the observed state variables
ObsX=[true];

% initial state variable for integration
X0=1;         
% initial filter estimates
Xf0=X0;       
% initial filter variance estimates
Pf0=0;        

% definition of the model equation 
ModEq = 'SE_ExactMul';
% definition of the observation equation
ObsEq = 'OE_ExactMul';
% definition of the exact predictions of the model
ExactSol = 'S_ExactMul';

FileNames=MakeFileNames(ModEq,ObsEq,ExactSol);

% absolute and relative tolerances for filtering
AbsTol=[1e-9,1e-9,1e-6];
RelTol=[1e-6,1e-6,1e-6];

