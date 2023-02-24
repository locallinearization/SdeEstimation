% Script for model definition: Non-autonomous linear SDE with two additive noises 
%      "Approximate linear minimum variance filters for continuous-discrete 
%      state space models: convergence and practical adaptive algorithms", 
%      Jimenez J.C., IMA Journal of Mathematical Control and Information, 
%      Volume 36, Issue 2, June 2019, Pages 341–378. 

deltat_min=0.0001;   % sampling period for integration
deltat=1;            % sampling period for filtering
t0=0.01;             % initial time 
tTotal=10+t0;        % sampling time  

% definition of state space parameters
Theta(1)=1/4; 
Theta(2)=5; 
Theta(3)=2; 
Theta(4)=0.1;  
Theta(5)=0.0001;     % observation noise variance

% index of the observed state variables
ObsX=[true];

% initial state variable for integration
X0=10;  
% initial filter estimates
Xf0=X0;      
% initial filter variance estimates
Pf0=0;         

% definition of the model equation 
ModEq = 'SE_ExactAdd';
% definition of the observation equation
ObsEq = 'OE_ExactAdd';
% definition of the exact predictions of the model
ExactSol = 'S_ExactAdd';

FileNames=MakeFileNames(ModEq,ObsEq,ExactSol);

% absolute and relative tolerances for filtering
AbsTol=[1e-9,1e-9,1e-6];
RelTol=[1e-6,1e-6,1e-6];



