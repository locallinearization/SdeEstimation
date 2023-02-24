function [Xf,Pf,Xp,Pp,Inn,PInn,mpts,LogLK]=LLfilter_Exact(t,Z,Xf0,Pf0,theta,FileNames,LKWindow)
%    Exact Minimum Variance Filter for the estimation of the unobserved componet x of the 
% continuous model
%        dx = f(t,x,theta) dt + g(t,x,theta) dw, 
% with dicrete observation equation 
%         Z = h(t,x,theta) + p(t,x,theta)n + e, 
% given a set of discrete observation of Z on t. 
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
%
% OutPuts:
%   Xf : LL filter  
%   Pf : LL filter variance 
%   Xp : LL prediction 
%   Pp : LL prediction variance  
%   Inn: Innovation
%   PInn: Innovation variance
%   mpts: number of misssing points used to compute the prediction between two observations (1 x n_m)
%   LogLK:  -2log(Likelihood)

%  References
%   - "Approximate linear minimum variance filters for continuous-discrete 
%      state space models: convergence and practical adaptive algorithms", 
%      Jimenez J.C.
%      IMA Journal of Mathematical Control and Information, 
%      Volume 36, Issue 2, June 2019, Pages 341�378. 

n_m=size(t,2);        % number of instants of times
n_v_e=size(Xf0,1);    % number of state variables
nobs=size(Z,1);       % number of observed variables

Xf=zeros(n_v_e,n_m);
Xp=zeros(n_v_e,n_m);
Pf=zeros(n_v_e,n_v_e,n_m);
Pp=zeros(n_v_e,n_v_e,n_m);
Inn=[zeros(nobs,1) inf*ones(nobs,n_m-1)];
PInn=zeros(nobs,nobs,n_m);

Zt=Xf0;
Vt=Pf0;
Xf(:,1)=Zt;
Xp(:,1)=Zt;
Pf(:,:,1) = Pf0;
Pp(:,:,1) = Pf0;

LogLK = 0.0;
[Prediction, Filter, ObsEq] = ModObsEqType(FileNames);
for k= 2:n_m
   
  Ztdt=Zt;
  Pt  =Vt + Zt*Zt';
  
  [Ztdt, Pt] = Prediction(t(k),t(k-1),Ztdt,Pt,theta);            % Exact predictions 
  Ut = Pt - Ztdt*Ztdt';
    
  [Zt, Vt, nu, Sv] = Filter(t(k),Ztdt,Ut,theta,Z(:,k),ObsEq);    % LL filter & innovation 
   
  if (LKWindow(1)<=t(k)) && (t(k)<=LKWindow(2))
     LogLK = LogLK + log(det(Sv)) + nu' /Sv * nu;                
  end
  
  Xp(:,k)   = Ztdt;    % LL prediction
  Pp(:,:,k) = Ut;      % LL prediction variance
  Xf(:,k)   = Zt;      % LL filter 
  Pf(:,:,k) = Vt;      % LL filter variance
  
  obs= find(Z(:,k)~=inf);
  Inn(obs,k)  = nu;      % innovation
  PInn(obs,obs,k) = Sv;  % innovation variance
   
end
mpts=zeros(1,n_m);

LogLK = LogLK + log(2.*pi).*n_m;



%----------------------------
% Defines the type of Model and Observation Equations

function [Prediction, Filter, ObsEq] = ModObsEqType(FileNames)

Prediction = FileNames{3};

ObsEq = FileNames{2};

switch nargout(ObsEq)
    case 3
        Filter=@Filter_Add;
    case 6
        Filter=@Filter_Mult;
    otherwise
        disp('Error in Observation definition')
end


