% Script for model definition: stochastic Fitzhugh-Nagumo model
% Example 1 in: 
%    "Estimation of distribution algorithms for the computation of innovation
%    estimators of diffusion processes" 
%    Z. Gonzalez, J.C. Jimenez, L. Lozada-Chang, R. Santana
%    Mathematics and Computers in Simulation 187 (2021) 449–467

deltat_min=0.0005;  % sampling period for integration
deltat=0.5;         % sampling period for filtering
t0=0;               % initial time  
tTotal=50;          % sampling time  

% definition of state space parameters
Theta(1)=1;
Theta(2)=1;
Theta(3)=0.1;
Theta(4)=0.001;   % observation noise variance for x1
Theta(5)=0.001;   % observation noise variance for x2

% index of the observed state variables
ObsX=[true; true];

% initial state variable for integration
X0=[-0.932309; -0.673218];
% initial filter estimates
Xf0=X0;                      
% initial filter variance estimates          
Pf0=0.01*ones(2,2);    

% definition of the model equation 
ModEq = 'SE_FitzhughNagumo';
% definition of the observation equation
ObsEq = 'OE_FitzhughNagumo';

FileNames=MakeFileNames(ModEq,ObsEq);

% absolute and relative tolerances for filtering
AbsTol=[1e-6,1e-6,0.001];
RelTol=[1e-3,1e-3,0.001];

