% Script for model definition: SIR Epidemic Model 
%    "A Contribution to the Mathematical Theory of Epidemics". Kermack W.O., McKendrick A.G., 
%    Proceedings of the Royal Society A: Mathematical, Physical and Engineering Sciences. 115 (1927) 700-721.

deltat=1;          % sampling period for filtering
t0=1;              % initial time
tTotal=t0+99;      % sampling time  

% definition of state space parameters
Theta(1) = 0.0004;
Theta(2) = 0.04;
Theta(3) = 1;           % observation noise variance of x2
Theta(4) = 1;           % observation noise variance of x3

% index of the observed state variables
ObsX=[false; true; true];

% initial state variable for integration;
X0=[997; 3; 0];   

% initial filter estimates
Xf0=X0;               
% initial filter variance estimates
Pf0=zeros(3);     

% definition of the model equation 
ModEq = 'SE_SIR';
% definition of the observation equation
ObsEq = 'OE_SIR';

FileNames=MakeFileNames(ModEq,ObsEq);

% absolute and relative tolerances for filtering
AbsTol=[1.0e-9,1.0e-9];
RelTol=[1.0e-6,1.0e-6];



