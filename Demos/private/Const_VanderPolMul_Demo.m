% Script for model definition:  Van der Pol model with random frecuency (multiplicative noise)
%     "Bias reduction in the estimation of diffusion processes from discrete observations",
%      Jimenez J.C, IMA Journal of Mathematical Control and Information, 
%      Volumen 37, December 2020, Pages 1468–1505. 

deltat_min=0.0001;   % sampling period for integration
deltat=1;            % sampling period for filtering
t0=0;                % initial time  
tTotal=t0+30;        % sampling time  

% definition of state space parameters
Theta(1)=1;
Theta(2)=1;          % frequency mean of the oscillator
Theta(3)=0;          % intensity mean of the random force
Theta(4)=1;          % system noise 
Theta(5)=0.001;      % observation noise variance

% index of the observed state variables
ObsX=[true; false];

% initial state variable for integration, X0(3) must be X0(1)^2
X0=[1; 1];          
% initial filter estimates
Xf0=X0;                 
% initial filter variance estimates
Pf0=zeros(2);             

% definition of the model equation 
ModEq = 'SE_VanderPolMul';
% definition of the observation equation
ObsEq = 'OE_VanderPolMul';

FileNames=MakeFileNames(ModEq,ObsEq);

% absolute and relative tolerances for filtering
AbsTol=[1e-6,1e-6,0.0001];
RelTol=[1e-4,1e-4,0.0001];
