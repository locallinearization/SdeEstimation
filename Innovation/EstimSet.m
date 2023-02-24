function EstimOptions = EstimSet(t,Z,Xf0,Pf0,theta0,FileNames,LKWindow,nu_stast,LB,UB,OptimMethod,OptimOptions,mpts,AbsTol,RelTol,NumSim)
%    EstimSet.m sets, for the parameter estimation, the information about the space state model with 
% continuous-time equation
%        dx = f(t,x,theta) dt + g(t,x,theta) dw, 
% and dicrete observation equation 
%         Z = h(t,x,theta) + p(t,x,theta) n + e, 
% given a set of discrete observation of Z. 
%
% Notation
%   n_m:   number of time-instants where Z is observed
%   n_obs: number of observation equations
%   n_v_e: number of state variables x 
%   n_theta: number of parameters theta 
%
% Input:
%   t: instants of times where at least one component of Z is observed  (1 x n_m)          
%   Z: observed sample                        (n_obs x n_m)   Set Z(k,t)=inf when data Z(k,t) is missing
%   Xf0: initial filter value                 (n_v_e x 1)
%   Pf0: initial filter variance value        (n_v_e x n_v_e)
%   theta0: initial parameter values for local optimization   (1 x n_theta)
%   FileNames={ModEq,ObsEq,ExactSol}   
%      ModEq string with the m file name of the continuous mnodel definition
%      ObsEq string with the m file name of the observation equation definition
%      ExactSol string with the m file name of the exact E(x/Zn) and E(xx'/Zn) 
%   LKWindow: interval of time where the likelihood function will be computed  (1 x 2)
%   nu_stast: a priori upper bound for the innovation mean and variance 
%   LB: lower bound for the parameters to be estimated (1 x n_v_e + n_theta)
%   UB: upper bound for the parameters to be estimated (1 x n_v_e + n_theta)
%   OptimMethod: name of the optimization method
%   OptimOptions: options of optimization method
%   mpts: (for non adaptive filter) number of a priori artificial points to compute the prediction between two observations (1 x n_m)
%   mpts: (for adaptive filter) if mpts=[], no restrinction to the maximum number of allowed missing points between the 
%                                          observations 
%                                  else,    mpts(k) indicates the maximum number of allowed missing points between the 
%                                          observations k-1 and k that could be adaptively estimated according 
%                                          to the tolerance especified in AbsTol(1:2) and RelTol(1:2)           (1 x n_m)    
%   AbsTol: Absolute Tolerance (1 x 3), 1: for the first moment, 
%                                       2: for the second moment, 
%                                       3: for the number of simulations
%   RelTol: Relative Tolerance (1 x 3), 1,2,3 same than AbsTol 
%   NumSim: Number of Simulations for the stochastic filter.
%           If NSim=[] not restriction on the maximum number of simulations
%           If NSim<>[], with AbsTol(3)=inf and RelTol(3)=inf, the stochastic adaptive LL filter performs NSim(k) 
%           simulations at each observation k.  (1 x n_m)
%
%   Set OptimOptions = [] to use dafaults optimization options
%   Set AbsTol = [] to use dafaults values AbsTol=[10^-6,10^-6,0.01] 
%   Set RelTol = [] to use dafaults values RelTol=[10^-3,10^-3,0.01]
%
%   For approximate innovation estimator with Deterministic Local Linearization Filters
%   [...] = EstimSet(t,Z,Xf0,Pf0,theta0,FileNames,LKWindow,nu_stast,LB,UB,OptimMethod,OptimOptions,mpts)
%         Sets the innovation estimator for using the non adaptive LL filter with the number of point between 
%         observations specified in input variable "mpts". 
%         For using the Classical LL filter, set mpts = zeros(1,n_m)
%   [...] = EstimSet(t,Z,Xf0,Pf0,theta0,FileNames,LKWindow,nu_stast,LB,UB,OptimMethod,OptimOptions,mpts,AbsTol,RelTol)
%         Sets innovation estimator for using the adaptive deterministic LL filter with the absolute 
%         and relative tolerances AbsTol,RelTol, and with a maximum of mpts(k) artificial points between 
%         the observation k-1 and k.
%         If AbsTol=RelTol=[], defaults AbsTol=[10^-6,10^-6], RelTol=[10^-3,10^-3] are used.
%         If mpts=[] not restriction on the number of artificial points between two observation
%
%   For approximate innovation estimator with Stochastic Local Linearization Filters
%   [...] = EstimSet(t,Z,Xf0,Pf0,theta0,FileNames,LKWindow,nu_stast,LB,UB,OptimMethod,OptimOptions,mpts,AbsTol,RelTol,NumSim)
%         Sets the innovation estimators for using the stochastic adaptive LL filter with the absolute 
%         and relative tolerances AbsTol,RelTol, with a maximum of mpts(k) artificial points between 
%         the observation k-1 and k, and a maximum of NSim(k) simulations at k.
%         If AbsTol=RelTol=[], defaults AbsTol=[10^-6,10^-6,0.01] and RelTol=[10^-3,10^-3,0.01] are used
%         If mpts=[] not restriction on the maximum number of artificial points between two observation
%         If NumSim=[] not restriction on the maximum number of simulations
%         If NumSim<>[], with AbsTol(3)=inf and RelTol(3)=inf, the stochastic adaptive LL filter performs NumSim(k) 
%            simulations at each observation k
%
%   For exact innovation estimator  
%   [...] = EstimSet(t,Z,Xf0,Pf0,theta0,FileNames,LKWindow,[],nu_stast,LB,UB,OptimMethod,OptimOptions,0,0)
%         Sets the innovation estimator for using the exact minimum variance filter with exact predictions between observations 

EstimOptions{1}=t;
EstimOptions{2}=Z;
EstimOptions{3}=Xf0;
EstimOptions{4}=Pf0;
EstimOptions{5}=theta0;
EstimOptions{6}=FileNames;
EstimOptions{7}=LKWindow;
EstimOptions{9}=nu_stast;
EstimOptions{10}=LB;
EstimOptions{11}=UB;
EstimOptions{12}=OptimOptions;
EstimOptions{16}=mpts;
EstimOptions{17}=OptimMethod;

switch nargin
   case 16 
     EstimOptions{20}=NumSim;
     EstimOptions{18}=AbsTol;
     EstimOptions{19}=RelTol;
     LLF = LLfilter(t,Z,Xf0,Pf0,theta0,FileNames,LKWindow,mpts,AbsTol,RelTol,NumSim);
   case 15
     EstimOptions{20}=-1;
     EstimOptions{18}=AbsTol;
     EstimOptions{19}=RelTol;
     LLF = LLfilter(t,Z,Xf0,Pf0,theta0,FileNames,LKWindow,mpts,AbsTol,RelTol);
    otherwise
     EstimOptions{20}=0;
     LLF = LLfilter(t,Z,Xf0,Pf0,theta0,FileNames,LKWindow,mpts);
end

if isnan(LLF.LogLK) || ~isreal(LLF.LogLK)
  error('Warning: Parameter initial values are unfeasible')
else    
  EstimOptions{8} = LLF.LogLK;   
  disp(['  '])
  % Log likelihood for the initial parameters Xf00 and theta0 (this is for normalization propose)
  disp([' Log LK for Theta0= ' num2str(-0.5*LLF.LogLK)])
end


