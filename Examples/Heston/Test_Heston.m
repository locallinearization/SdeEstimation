% Script for state estimation
close all

% Model definition
Const_Heston;

load data_Heston
Z=datos;

% plot state variables and observations
PlotModelRealization(t,Z,proceso,ObsX);

% interval of time where the likelihood function will be computed
LKWindow=[t(1) t(end)];

% setting the LL filter
% FilterMethod='ClassFilter';         % Classical LL filter
% FilterMethod='DeterFilter';         % Deterministic LL filter and a priori artificial points between observations
% FilterMethod='StochFilter';         % Stochastic LL Filter
  FilterMethod='AdapDeterFilter';     % Adaptive Deterministic LL Filter
switch FilterMethod
    case 'ClassFilter'
       n_m = length(t);
       mpts = zeros(1,n_m);       % classical LL filter
       LLF = LLfilter(t,Z,Xf0,Pf0,Theta,FileNames,LKWindow,mpts);
    case 'DeterFilter'
       n_m = length(t);
       mpts = 3.*ones(1,n_m);     % a priori artificial points between observation times
       LLF = LLfilter(t,Z,Xf0,Pf0,Theta,FileNames,LKWindow,mpts);
    case 'StochFilter'
       n_m = length(t);
       mpts = 3.*ones(1,n_m);   % a priori artificial points between observation times
       AbsTol=[0,0,1e-7];
       RelTol=[0,0,1e-7];
       NumSim=[];
       LLF = LLfilter(t,Z,Xf0,Pf0,Theta,FileNames,LKWindow,mpts,AbsTol,RelTol,NumSim);      
    case 'AdapDeterFilter'
       mpts=[];
       LLF = LLfilter(t,Z,Xf0,Pf0,Theta,FileNames,LKWindow,mpts,AbsTol,RelTol);
    otherwise
       error('Filter Method is not well defined');    
end

disp(['   '])
disp([' Log LK = ' num2str(-LLF.LogLK/2)])

% plot observations and filter estimates
PlotFilterEstimates(t,Z,LLF,ObsX)
