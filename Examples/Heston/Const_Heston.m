% Script for model definition: Heston model
%   "A closed-form solution for options with stochastic volatility with applications
%    to bond and currency options", Heston, S. L. (1993) , Rev. Fin. Stud. 6, 327–343. 
% Model setting from: 
%   "State and parameter estimation of stochastic physical systems 
%    from uncertain and indirect measurements", Jimenez J.C., Yoshimoto A.,
%    Miwakeich F., The European Physical Journal Plus, 136 (2021) 1-17. 

deltat_min=1/25000;   % sampling period for integration
deltat=1/250;         % sampling period for filtering (daily sampling)
t0=0;                 % initial time  
tTotal=t0+17;         % sampling time  

% Heston model parameters
alpha = 0.067;        % expected return
lamda = 1.980;        % rate at which the variance reverts to its long-term mean,
mu = 0.0197;          % mean long-term variance
gamma = 0.069;        % the volatility of the variance process.
ro =  -0.491;         % correlation between the Wiener processes

% definition of state space parameters for the transformed Heston model
Theta(1)= alpha;       
Theta(2)= mu*lamda;    
Theta(3)= lamda;     
Theta(4)= gamma*ro;     
Theta(5)= gamma*sqrt(1-ro^2);        
Theta(6)=0.0001;     % observation noise variance

% index of the observed state variables
ObsX=[true; false];

% initial state variable for integration
X0=[log(353); 0.01];          
% initial filter estimates
Xf0=X0;                 
% initial filter variance estimates
Pf0=diag([0.001 0.000001]);             

% definition of the model equation 
ModEq = 'SE_Heston';
% definition of the observation equation
ObsEq = 'OE_Heston';

FileNames=MakeFileNames(ModEq,ObsEq);

% absolute and relative tolerances for filtering
AbsTol=[1e-6,1e-6,0.001];
RelTol=[1e-3,1e-3,0.001];



