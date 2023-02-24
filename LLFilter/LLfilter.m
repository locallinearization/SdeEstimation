function varargout = LLfilter(t,Z,Xf0,Pf0,theta,FileNames,LKWindow,mpts,AbsTol,RelTol,NumSim)
%    Local Linearization Filters for the estimation of the unobserved componet x of the 
% continuous model
%        dx = f(t,x,theta) dt + g(t,x,theta) dw, 
% with dicrete observation equation 
%         Z = h(t,x,theta) + p(t,x,theta) n + e, 
% given a set of discrete observation of Z on t. 
%
% Xf0, Pf0 and theta are also known and given in advance.
%
% Notation
%   n_m: number of time-instants where Z is observed
%   n_obs: number of observation equations
%   n_v_e: number of state variables x 
%   n_theta: number of parameters theta 
%
% Input:
%   t: instants of times where at least one component of Z is observed  (1 x n_m)          
%   Z: observed sample                        (n_obs x n_m)   Set Z(k,t)=inf when data Z(k,t) is missing
%   Xf0: initial filter value                 (n_v_e x 1)
%   Pf0: initial filter variance value        (n_v_e x n_v_e)
%   theta: State Space Model parameters       (1 x n_theta)
%   FileNames={ModEq,ObsEq,ExactSol}   
%      ModEq string with the .m file name of the continuous model definition
%      ObsEq string with the .m file name of the observation equation definition
%      ExactSol string with the .m file name of the exact E(x/Zn) and E(xx'/Zn) 
%   LKWindow: interval of time where the likelihood function will be computed  (1 x 2)
%             set LKWindow=[] to get LKWindow = [t(1) t(end)]
%   mpts: (for non adaptive filter) number of a priori artificial points to compute the prediction between two observations (1 x n_m)
%   mpts: (for adaptive filter) if mpts=[], no restrinction to the maximum number of allowed artificial points between the 
%                                          observations 
%                                  else,    mpts(k) indicates the maximum number of allowed artificial points between the 
%                                          observations k-1 and k that could be adaptively estimated according 
%                                          to the tolerance especified in AbsTol(1:2) and RelTol(1:2)           (1 x n_m)    
%   AbsTol: Absolute Tolerance (1 x 3), 1: for the first moment, 
%                                       2: for the second moment, 
%                                       3: for the number of simulations
%   RelTol: Relative Tolerance (1 x 3), 1,2,3 same than AbsTol 
%   NumSim: Maximun Number of Simulations for the stochastic filter.
%           If NSim=[] not restriction on the maximum number of simulations
%           If NSim<>[], with AbsTol(3)=inf and RelTol(3)=inf, the stochastic adaptive LL filter performs NSim(k) 
%           simulations at each observation k.  (1 x n_m)
%
% OutPuts:
%   Xf : LL filter  
%   Pf : LL filter variance 
%   Xp : LL prediction 
%   Pp : LL prediction variance  
%   Inn: Innovation
%   PInn: Innovation variance
%   pts: (deterministic filter) number of artificial points used to compute 
%         the prediction between two observations (2 x n_m)
%         row 1: accepted, row 2: failed
%   pts: (stochastic filter) average (per simulations) of artificial points used to compute 
%         the prediction between two observations (3 x n_m)
%         row 1: accepted, row 2: failed, row 3: total number of simulations  
%   LogLK:  -2log(Likelihood)
%
% Call options for Input variables:
%
%   Deterministic Local Linearization Filters
%     [...] =  LLfilter(t,Z,Xf0,Pf0,theta,FileNames,LKWindow) computes the 
%     adaptive LL filter with the default absolute and relative tolerance 10^-6 and 10^-3
%     [...] =  LLfilter(t,Z,Xf0,Pf0,theta,FileNames,LKWindow,mpts) computes the 
%     non adaptive LL filter with number of points between observations specified in the input
%     variable "mpts". For Classical LL filter, set mpts = zeros(1,n_m)
%     [...] =  LLfilter(t,Z,Xf0,Pf0,theta,FileNames,LKWindow,mpts,AbsTol,RelTol) computes the
%     adaptive LL filter with the specified absolute and relative tolerances AbsTol,RelTol, 
%     with a maximum of mpts(k) artificial points between the observation k-1 and k.
%     If AbsTol=RelTol=[], defaults AbsTol=[10^-6,10^-6], RelTol=[10^-3,10^-3] are used.
%     If mpts=[] not restriction on the number of artificial points between two observation
%                
%   Stochastic Local Linearization Filters
%     [...] =  LLfilter(t,Z,Xf0,Pf0,theta,FileNames,LKWindow,mpts,AbsTol,RelTol,NumSim), computes the
%     stochastic LL filter with mpts(k) artificial points between the observation k-1 and k, 
%     and a maximum of NumSim(k) simulations at k.
%     If NumSim=[] not restriction on the maximum number of simulations
%     If NumSim<>[], with AbsTol(3)=inf and RelTol(3)=inf, the stochastic LL filter performs 
%        exactly NumSim(k) simulations at each observation k
%     
%   Exact Minimum Variance Filter
%      [...] =  LLfilter(t,Z,Xf0,Pf0,theta,FileNames,LKWindow,[],0,0) computes the 
%      exact minimum variance filter with exact predictions between observations 
%
% Call options for Output variables:
%
%   for more than one output
%      [Xf,Pf,Xp,Pp,Inn,PInn,pts,LogLK] =  LLfilter(...)
%
%   for one output
%      LLF =  LLfilter(...), where LLF is a record variable with the following structure
%                LLF.Xf 
%                LLF.Pf
%                LLF.Xp
%                LLF.Pp
%                LLF.Inn
%                LLF.PInn
%                LLF.pts
%                LLF.LogLK
%                LLF.AlertPts: 1 if the maximun number of artificial point between observation is reached 
%                LLF.AlertSim: 1 if the maximun number of simulations is reached

%  References
%   - "Approximate linear minimum variance filters for continuous-discrete 
%      state space models: convergence and practical adaptive algorithms", 
%      Jimenez J.C.
%      IMA Journal of Mathematical Control and Information, 
%      Volume 36, Issue 2, June 2019, Pages 341–378. 
%   - "Simplified formulas for the mean and variance of linear 
%      stochastic differential equations",  
%      Jimenez J.C.
%      Applied Mathematics Letters
%      Volume 49, November 2015, Pages 12-19. 

if nargin==7
   Atol=[1.0e-6,1.0e-6];
   Rtol=[1.0e-3,1.0e-3];
   mpts=[];
end

if (nargin>=9) && (~isempty(AbsTol)) 
   if length(AbsTol)==2 || length(AbsTol)==3
     Atol=AbsTol;     
   else
     Atol=[AbsTol,AbsTol];     
   end
else
   Atol=[1.0e-6,1.0e-6];
end

if (nargin>=10) && (~isempty(RelTol)) 
   if length(RelTol)==2 || length(RelTol)==3
     Rtol=RelTol; 
   else
     Rtol=[RelTol,RelTol];     
   end
else
   Rtol=[1.0e-3,1.0e-3];
end

if (nargin>=11) && length(RelTol)==2
      Atol=[Atol 0.01];
      Rtol=[Rtol 0.01];
end

if isempty(LKWindow)
    LKWindow=[t(1) t(end)];
end


switch nargin
    case 11
      % Non-Adaptive time-stepping Stochastic Local Linearization Filter
      [Xf,Pf,Xp,Pp,Inn,PInn,pts,LogLK]=LLfilter_SimMpts(t,Z,Xf0,Pf0,theta,FileNames,LKWindow,mpts,Atol,Rtol,NumSim);
    case {7,10}
       if (sum(Atol)==0) && (sum(Rtol)==0) && isempty(mpts)       
          % Exact Minimum Variance Filter
          [Xf,Pf,Xp,Pp,Inn,PInn,pts,LogLK]=LLfilter_Exact(t,Z,Xf0,Pf0,theta,FileNames,LKWindow);
       else 
          % Adaptive time-stepping Deterministic Local Linearization Filter
          [Xf,Pf,Xp,Pp,Inn,PInn,pts,LogLK]=LLfilter_Det(t,Z,Xf0,Pf0,theta,FileNames,LKWindow,mpts,Atol,Rtol);
       end
    case 8
       % Non-Adaptive time-stepping Deterministic Local Linearization Filter
       [Xf,Pf,Xp,Pp,Inn,PInn,pts,LogLK]=LLfilter_DetMpts(t,Z,Xf0,Pf0,theta,FileNames,LKWindow,mpts);
end

if nargout == 1
  varargout{1}.Xf=Xf;
  varargout{1}.Pf=Pf;
  varargout{1}.Xp=Xp;
  varargout{1}.Pp=Pp;
  varargout{1}.Inn=Inn;
  varargout{1}.PInn=PInn;
  varargout{1}.pts=pts;
  varargout{1}.LogLK=LogLK;
  if ~((sum(Atol(1:2))==0) && (sum(Rtol(1:2))==0) && ~isempty(mpts))       
     varargout{1}.AlertPts = ~isempty(mpts) && (size(pts,1)>1) && (sum(sum(pts(1:2,:))>=mpts)>0);  
  end
  if (nargin>=11)
    if (length(NumSim)>1) && (AbsTol(3)==inf) && (RelTol(3)==inf)
      varargout{1}.AlertSim = 0;
    else
      varargout{1}.AlertSim = (length(NumSim)>1) && (sum(sum(pts(3,:))>=NumSim)>0);  
    end 
  end
else
  varargout{1}=Xf;
  varargout{2}=Pf;
  varargout{3}=Xp;
  varargout{4}=Pp;
  varargout{5}=Inn;
  varargout{6}=PInn;
  varargout{7}=pts;
  varargout{8}=LogLK;
end
