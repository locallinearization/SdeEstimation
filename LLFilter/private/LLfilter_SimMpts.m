function [Xf,Pf,Xp,Pp,Inn,PInn,pts,LogLK]=LLfilter_SimMpts(t,Z,Xf0,Pf0,theta,FileNames,LKWindow,mpts,AbsTol,RelTol,NumSim)
%    Non-Adaptive Stochastic Local Linearization Filter for the estimation of the unobserved componet x of the 
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
%   mpts: number of misssing points to compute the prediction between two observations (1 x n_m)
%   AbsTol: Absolute Tolerance (1 x 3), 1: for the first moment, 
%                                       2: for the second moment, 
%                                       3: for the number of simulations
%   RelTol: Relative Tolerance (1 x 3), 1,2,3 same than AbsTol 
%   NumSim: Number of Simulations
%           If NumSim=[], NumSim is adaptively estimated according to the tolerance especified in AbsTol(3) and RelTol(3)
%           If NumSim<>[], with AbsTol(3)=inf and RelTol(3)=inf, the stochastic adaptive LL filter performs NumSim(k) 
%           simulations at each observation k
%
% OutPuts:
%   Xf : LL filter  
%   Pf : LL filter variance 
%   Xp : LL prediction 
%   Pp : LL prediction variance  
%   Inn: Innovation
%   PInn: Innovation variance
%   pts: (stochastic filter) number of misssing points used to compute 
%         the prediction between two observations (3 x n_m)
%         row 1: mpts, row 2: zeros, row 3: total number of simulations  
%   LogLK:  -2log(Likelihood)
%

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

if ~isempty(mpts) && (length(mpts)~=length(t))
    error('max_pts wrong. Check the lenght of max_pts')
end
if ~isempty(NumSim) && (length(NumSim)~=length(t))
    error('NumSim wrong. Check the lenght of NumSim')
end

n_m=size(t,2);        % number of instants of times
n_v_e=size(Xf0,1);    % number of state variables
nobs=size(Z,1);       % number of observed variables

LocalCluster=parcluster;
NW=LocalCluster.NumWorkers;
MinSim=NW;
while MinSim<100
    MinSim=2*MinSim;
end
if ~isempty(NumSim)
    MinSim=min([MinSim min(NumSim)]);
end

AdaptiveSampling = (RelTol(3)~=inf) && (AbsTol(3)~=inf); 
if AdaptiveSampling
  TNSim=MinSim*ones(1,n_m); 
  if isempty(NumSim), NumSim=inf*ones(1,n_m); end 
else
  TNSim=NumSim;
end
TNSim(1)=1;

Xf=zeros(n_v_e,n_m);
Xp=zeros(n_v_e,n_m);
Pf=zeros(n_v_e,n_v_e,n_m);
Pp=zeros(n_v_e,n_v_e,n_m);
Inn=[zeros(nobs,1) inf*ones(nobs,n_m-1)];
PInn=zeros(nobs,nobs,n_m);
uno = ones(n_v_e,1);

Zt=Xf0;
Vt=Pf0;
Xf(:,1)=Zt;
Xp(:,1)=Zt;
Pf(:,:,1) = Pf0;
Pp(:,:,1) = Pf0;

LogLK = 0.0;
[Prediction, Filter, ModEq, ObsEq] = ModObsEqType(FileNames,n_v_e);

% this is for controling the Monte Carlo simulation
SemillaMC = 1000*reshape(1:NW*n_m,NW,n_m);

if AdaptiveSampling
  GCF=figure(10000);
  XPos=GCF.Position(1)-390;
  YPos=GCF.Position(2)+GCF.Position(4)-200;
  close(GCF)
  fig = uifigure('Name','Stochastic simulation','Position',[XPos YPos 380 250]);
  uilabel(fig,'Position',[10 fig.Position(4)-20 500 15],'Text','Obs   StatErrorMean   StatErrorVar   SimMean   SimVar   Sim');
  txa = uitextarea(fig,'Position',[10 10 fig.Position(3)-20 fig.Position(4)-35]);
end
for k= 2:n_m
  Y=[];
  NSim=TNSim(k);
  NoEnd=1; 
  while NoEnd 
    NSimPool=ones(1,NW)*round(NSim/NW);
    NSimPool(end)=NSimPool(end)+mod(NSim,NW);

    spmd
        rand('seed',SemillaMC(labindex,k));   % this is for controling the Monte Carlo simulation
        
        y=zeros(NSimPool(labindex),n_v_e);
        for n=1:NSimPool(labindex)
                    
           [U,D,V] = svd(Vt);
           S = U*sqrt(D)*V';
           Ruido = sign(2*rand(n_v_e,1)-uno);
           Ztdt = Zt + S*Ruido;
           Pt   = Ztdt*Ztdt';
                    
           hh  =(t(k)-t(k-1))/(mpts(k)+1);
           for j=1:(mpts(k)+1)                                                   % LL prediction
              [Ztdt, Pt] = Prediction(t(k-1)+j*hh,t(k-1)+(j-1)*hh,Ztdt,Pt,theta,ModEq,1.0e-9,1.0e-12);  
              
              W = Pt - Ztdt*Ztdt';
              [U,D,V] = svd(W);
              S = U*sqrt(D)*V';
              Ruido = sign(2*rand(n_v_e,1)-uno);

              Ztdt = Ztdt + S*Ruido;
              Pt   = Ztdt*Ztdt';
           end

           y(n,:) = Ztdt;
        end
    end
    for p=1:NW
        Y=[Y; y{p}];
    end
    if AdaptiveSampling 
       if NumSim(k)==TNSim(k)
           [TNSim(k), NoEnd, Info]=num_samples(M1,M2,M4,TNSim(k),RelTol(3),AbsTol(3),NW,NumSim(k));
           if isvalid(fig)
             txa.Position(3:4)=[fig.Position(3)-20 fig.Position(4)-30]; 
             txa.Value=[txa.Value(1:k-1); {[num2str(k),'     ',Info]}];
           end
       else         
         M1=mean(Y)';
         M2=mean(Y.^2)';
         M4=mean(Y.^4)';
         NSim=TNSim(k);
         [TNSim(k), NoEnd, Info]=num_samples(M1,M2,M4,TNSim(k),RelTol(3),AbsTol(3),NW,NumSim(k));
         if isvalid(fig)
           txa.Position(3:4)=[fig.Position(3)-20 fig.Position(4)-30];
           if NoEnd
             txa.Value=[txa.Value; {[num2str(k),'     ',Info]}];
           else
             txa.Value=[txa.Value(1:k-1); {[num2str(k),'     ',Info]}];
           end
         end
         NSim=TNSim(k)-NSim;
       end
       drawnow
    else
      NoEnd=0;
    end
  end
  Ztdt=mean(Y)';
  Ut = cov(Y);
 
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
  PInn(obs,obs,k) = Sv;      % innovation variance
   
end

pts=[mpts; zeros(1,n_m); TNSim];

LogLK = LogLK + log(2.*pi).*n_m;

if AdaptiveSampling, delete(fig), end

