% Script for computing the Innovation Estimate of the unknown parameters "theta" and the unobserved componet x of the 
% continuous-time model
%        dx = f(t,x,theta) dt + g(t,x,theta) dw, 
% with dicrete observation equation 
%         Z = h(t,x,theta) + p(t,x,theta) n + e, 
% given a set of discrete observation of Z. The initial filter value Xf0 can be also estimated.
%
% For Z = x, the Innovation estimator of "theta" reduces to the Quasi Maximum Likelihood estimator of "theta"
% 
% Four optimization algorithms are available: 
% - three locals ('fmincon','fminsearch','patternsearch') 
% - a global ('UMDA+fmincon')
%
% Five options of filtering algorithms are available: 
% - Local Linearization Filters (classical and with prior artificial points between observation times)
% - Adaptive Deterministic Local Linearization Filters 
% - Adaptive Stochastic Local Linearization Filters 
% - Exact Minimum Variance Filter, when the exact predictions are known

close all

% Model definition
Const_Heston;

load data_Heston
Z=datos;

%*********************  Estimation of the unknown model parameters (theta)

% initial state variable to be estimated
LB_X=[];        % lower bound
UB_X=[];        % upper bound
Xf0Index=[];

% parameters to be estimated
LB_T = [0   0.01 1 -0.2  0    0.00001];     
UB_T = [0.2 0.05 2  0    0.3  0.0005];  
ThetaIndex=[1 2 3 4 5];

% state initial values for local optimization
Xf0=X0;

% parameter initial values for local optimization
Theta0=Theta;     

% lower and upper bound of the parameters to be estimated
LB=[LB_X(Xf0Index) LB_T(ThetaIndex)];
UB=[UB_X(Xf0Index) UB_T(ThetaIndex)];

% setting the optimization method
OptimMethod='fmincon';
% OptimMethod='fminsearch';
% OptimMethod='patternsearch';
% OptimMethod='UMDA+fmincon';
switch OptimMethod
  case 'fmincon'
    TolFunOpt=1e-6;                
    TolXOpt=1e-6;            
    MaxFunEval = 300*length(ThetaIndex);
    OptimOptions = optimset('fmincon');
    OptimOptions = optimset(OptimOptions,'LargeScale','off','Display','iter','MaxFunEval',1e3,'TolFun',TolFunOpt,'TolX',TolXOpt);
  case 'fminsearch'
    TolFunOpt=1e-6;                
    TolXOpt=1e-6;            
    OptimOptions = optimset('fminsearch');
    OptimOptions = optimset(OptimOptions,'Display','iter','MaxFunEval',1e3,'TolFun',TolFunOpt,'TolX',TolXOpt);
  case 'patternsearch'
    TolFunOpt=1e-6;                
    TolXOpt=1e-3;            
    OptimOptions = psoptimset;
    %OptimOptions = psoptimset(OptimOptions,'PollMethod','MADSPositiveBasis2N','Display','iter','TolFun',TolFunOpt,'TolX',TolFunOpt);
    OptimOptions = psoptimset(OptimOptions,'Display','iter','TolFun',TolFunOpt,'TolX',TolFunOpt);
  case 'UMDA+fmincon'
    % definitions for fmincom
    TolFunOpt=1e-6;                
    TolXOpt=1e-6;            
    OptimOptions = optimset('fmincon');
    OptimOptions = optimset(OptimOptions,'LargeScale','off','Display','iter','MaxFunEval',1e3,'TolFun',TolFunOpt,'TolX',TolXOpt);
    % definitions for UMDA
    OptimOptions.PopSizeRate = 20;    % PopSizeRate * length(ParVar) = Population Size;
    OptimOptions.NumGen = 10;         % Number of Generations
    OptimOptions.TruncValue = 0.5;    % Truncation value
    OptimOptions.MaxFitness = inf;    % Maximum value of the fitness function, useful as stop condition
  otherwise
    error('Optimization Method is not well defined');    
end

% optimization with a priory statistics on the innovations
nu_stast=[];

% interval of time where the likelihood function will be computed
LKWindow=[t(1) t(end)];

% setting the innovation estimator
% InnovMethod='ClassInnov';         % with Classical LL filter
% InnovMethod='DeterInnov';         % with Deterministic LL filter and a priori artificial points between observations
% InnovMethod='StochInnov';         % with Stochastic LL Filter, but with a priori artificial points between observations
  InnovMethod='AdapDeterInnov';     % with Adaptive Deterministic LL Filter
switch InnovMethod
    case 'ClassInnov'
       n_m = length(t);
       mpts = zeros(1,n_m);       % classical LL filter
       EstimOptions = EstimSet(t,Z,Xf0,Pf0,Theta0,FileNames,LKWindow,nu_stast,LB,UB,OptimMethod,OptimOptions,mpts);
    case 'DeterInnov'
       n_m = length(t);
       mpts = 10.*ones(1,n_m);    % a priori artificial points between observation times
       EstimOptions = EstimSet(t,Z,Xf0,Pf0,Theta0,FileNames,LKWindow,nu_stast,LB,UB,OptimMethod,OptimOptions,mpts);
    case 'StochInnov'
       n_m = length(t);
       mpts = 3.*ones(1,n_m);   % a priori artificial points between observation times
       NumSim=[];
       AbsTol=[0,0,1e-7];
       RelTol=[0,0,1e-7];
       EstimOptions = EstimSet(t,Z,Xf0,Pf0,Theta0,FileNames,LKWindow,nu_stast,LB,UB,OptimMethod,OptimOptions,mpts,AbsTol,RelTol,NumSim);
    case 'AdapDeterInnov'
       mpts=[];
       EstimOptions = EstimSet(t,Z,Xf0,Pf0,Theta0,FileNames,LKWindow,nu_stast,LB,UB,OptimMethod,OptimOptions,mpts,AbsTol,RelTol);
    otherwise
       error('Innovation Method is not well defined');    
end

% computing the innovation estimators
InnEst = InnovEstimator(Xf0Index,ThetaIndex,EstimOptions);

% Summary of the parameter estimation
DisplayParameterEstimates(InnEst,X0,Theta,Xf0Index,ThetaIndex,EstimOptions)

% plot observations and fitting-filter estimates
PlotFilterEstimates(t,Z,InnEst.LLF,ObsX)

% plot statistics of the standarized fitting-innovation
PlotStatTestFitInn(InnEst.LLF)

% Recover original Heston parameter values
RecoverHestonParameters


