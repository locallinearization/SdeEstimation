% Script for model definition: Harmonic Oscillator with random frecuency and force
%   "Simplified formulas for the mean and variance of linear stochastic 
%    differential equations", Jimenez J.C., Appl. Math. Letters, 49 (2015) 12-19  
% Model setting from: 
%   "State and parameter estimation of stochastic physical systems 
%    from uncertain and indirect measurements", Jimenez J.C., Yoshimoto A.,
%    Miwakeich F., The European Physical Journal Plus, 136 (2021) 1-17.

deltat_min=0.0001;   % sampling period for integration
deltat=0.1;          % sampling period for filtering
t0=0;                % initial time  
tTotal=t0+80;        % sampling time  

% definition of state space parameters
Theta(1)=10;        % intensity mean of the random force
Theta(2)=30;        % frequency mean of the oscillator
Theta(3)=0.5;       % noisy frecuency
Theta(4)=0.5;       % noisy force
Theta(5)=0.1;       % observation noise variance

% index of the observed state variables
ObsX=[false; true];

% initial state variable for integration
X0=[1; 1];          
% initial filter estimates
Xf0=X0;                 
% initial filter variance estimates
Pf0=zeros(2);

% definition of the model equation 
ModEq = 'SE_LinearOscillator';
% definition of the observation equation
ObsEq = 'OE_LinearOscillator';
% definition of the exact predictions of the model
ExactSol = 'S_LinearOscillator';

FileNames=MakeFileNames(ModEq,ObsEq,ExactSol);

% absolute and relative tolerances for filtering
AbsTol=[1e-6,1e-6,0.001];
RelTol=[1e-3,1e-3,0.001];

