function varargout = InnovEstimator(Xf0Index,ThetaIndex,EstimOptions)
%    Innovation estimator of the unknown parameters theta and unobserved componets x of the 
% continuous-time model
%        dx = f(t,x,theta) dt + g(t,x,theta) dw, 
% with dicrete observation equation 
%         Z = h(t,x,theta) + p(t,x,theta) n + e, 
% given a set of discrete observation of Z. The initial filter value Xf0 can be also estimated.
%
% For Z = x, the Innovation estimator of theta reduces to the Quasi Maximum Likelihood estimator of theta
%
% Notation
%   n_m:   number of time-instants where Z is observed
%
% Variable InPuts
%   Xf0Index:     Indexes of the Xf0 to be estimates ([] is no one)
%   ThetaIndex:   Indexes of the theta to be estimates ([] is no one)
%   EstimOptions: Full infomation about the Space State Model to be estimated (see EstimSet.m)
%
% Command Window OutPut during the estimation
%   Current state of the optimization algorithm
%
% Grafical OutPuts during the estimation
%   Standarized fitting-innovation inside the optimization algorithm
%   Number of accepted artificial points (apts) between two consecutive observations
%   Number of failed artificial points (fpts) between two consecutive observations
%   Number of stochastic simulations (NSim) between two consecutive observations
%
% Additional Window OutPut during the estimation with stochastic filter
%   Current Statistical Error (StatErrorMean and StatErrorVar) in the computation 
%   of the prediction mean and variance at each observation (Obs). Additional Number of
%   Simulations (SimMean and SimVar) needed to compute the prediction mean and variance
%   with the required tolerance. Total number of performed simulations(Sim).
%
% Variable OutPuts:
%   Xf0:         innovation estimates of Xf0
%   Theta:       innovation estimates of theta
%   LLF:         LL fitting filter estimates
%   Xf0CI:       90% confidential interval of Xf0
%   ThetaCI:     90% confidential interval of Theta
%   Xf0Cov:      Covariance matrix of Xf0
%   ThetaCov:    Covariance matrix of Theta
%   AIC:         Akaike Information Criterion
%   BIC:         Bayesian Information Criterion
%   OptExitFlag: exit condition of the optimizacion process
%   OptOuPut:    information about the end of the optimization process
% 
%   Here, LLF is a record variable with the following structure
%      LLF.Xf:       LL fitting-filter
%      LLF.Pf:       LL fitting-filter variance
%      LLF.Xp:       LL fitting-prediction
%      LLF.Pp:       LL fitting-prediction variance
%      LLF.Inn:      LL fitting-innovation
%      LLF.PInn:     LL fitting-innovation variance
%      LLF.pts:     (deterministic filter) number of a priori artificial points used to compute 
%                    the fitting-predictions between two observations (2 x n_m)
%                    row 1: accepted, row 2: failed
%      LLF.pts:     (stochastic filter) number of a priori artificial points used to compute 
%                    the fitting-predictions between two observations (3 x n_m)
%                    row 1: accepted, row 2: failed, row 3: total number of simulations  
%      LLF.LogLK:    -2log(Likelihood)
%      LLF.AlertPts: 1 if the maximun number of artificial point between observation is reached. 0, otherwise
%      LLF.AlertSim: 1 if the maximun number of simulations is reached. 0, otherwise
%      LLF.KStest:   Kolmogorov-Smirnov test of normality N(0,1) for the standarized fitted-innovation (record variable)
%      LLF.JBtest:   Jarque-Bera test of composite Gaussianity for the standarized fitted-innovation (record variable)
%      LLF.LBQtest:  Ljung-Box Q-test for the autocorrelation of the standarized fitted-innovation (record variable)
%      LLF.ARCHtest: Eagle ARCH test of heteroscedasticity for the standarized fitted-innovation (record variable)
%
%  JBtest Structure
%         JBtest.h      - Boolean decisions for the test. h equal to 1 indicates rejection of the
%                         null of Gaussianity. h equal to 0 indicates a failure to reject the null.
%         JBtest.pValue - p-values of the test statistic
%         JBtest.stat   - test statistic
%         JBtest.cValue - critical values for the tests at 90%
%
%  KStest Structure
%         Similar to the JBtest structure
%
%  LQBtest Structure
%         LQBtest.h -      Vector of Boolean decisions for the tests computed at the lags 0,...,min([20 n_m]) 
%   	                   Values of h equal to 1 indicate rejection of the null of hypotesis in favor
%                          of the alternative. Values of h equal to 0 indicate a failure to reject the null.
%         LQBtest.pValue - Vector of p-values of the test statistics at the mentioned lags   
%         LQBtest.stat -   Vector of test statistics at the mentioned lags
%         LQBtest.cValue - Vector of 90% critical values for the tests at the mentioned lags
%
%  ARCHtest Structure
%         Similar to the LQBtest structure
%
%  OptOutPut Structute
%         It is provided by the optimization algorithm used by the
%         innovation method
%
% A call to the function estimset.m MUST be done before to use InnovEstimators.m
%
% Call options for Input variables:
%
%   Approximate Innovation Estimator with Classical (deterministic) Local Linearization Filter
%      mpts = zeros(1,n_m);    
%      EstimOptions = EstimSet(t,Z,Xf0,Pf0,Theta0,FileNames,LKWindow,nu_stast,LB,UB,OptimMethod,OptimOptions,mpts);
%      [...] = InnovEstimators(Xf0Index,ThetaIndex,EstimOptions)
%
%   Approximate Innovation Estimator with non adaptive Deterministic Local Linearization Filter
%      mpts = a vector of n_m elements specifying the number of artificial points between consecutive observations 
%      EstimOptions = EstimSet(t,Z,Xf0,Pf0,Theta0,FileNames,LKWindow,nu_stast,LB,UB,OptimMethod,OptimOptions,mpts);
%      [...] = InnovEstimators(Xf0Index,ThetaIndex,EstimOptions)
%
%   Approximate Innovation Estimator with adaptive Deterministic Local Linearization Filter (filter's options in LLfilter.m apply here)
%      EstimOptions = EstimSet(t,Z,Xf0,Pf0,Theta0,FileNames,LKWindow,nu_stast,LB,UB,OptimMethod,OptimOptions,mpts,AbsTol,RelTol);
%      [...] = InnovEstimators(Xf0Index,ThetaIndex,EstimOptions)
%
%   Approximate Innovation Estimator with Stochastic Local Linearization Filter (filter's options in LLfilter.m apply here)
%      EstimOptions = EstimSet(t,Z,Xf0,Pf0,Theta0,FileNames,LKWindow,nu_stast,LB,UB,OptimMethod,OptimOptions,mpts,AbsTol,RelTol,NumSim);
%      [...] = InnovEstimators(Xf0Index,ThetaIndex,EstimOptions)
%
%   Exact Innovation Estimator  
%      EstimOptions = EstimSet(t,Z,Xf0,Pf0,Theta0,FileNames,LKWindow,mpts,nu_stast,LB,UB,OptimMethod,OptimOptions,[],0,0);
%      [...] = InnovEstimators(Xf0Index,ThetaIndex,EstimOptions)
%
% Call options for Output variables:
%
%   for more than one output
%      [Xf0,Theta,LLF,Xf0CI,ThetaCI,Xf0Cov,ThetaCov,AIC,BIC] =  InnovEstimators(...)
%
%   for one output
%      InnEst =  InnovEstimators(...), where InnEst is a record variable with the following structure
%                InnEst.Xf0     
%                InnEst.Theta   
%                InnEst.LLF     
%                InnEst.Xf0CI    
%                InnEst.ThetaCI  
%                InnEst.Xf0Cov   
%                InnEst.ThetaCov 
%                InnEst.AIC
%                InnEst.BIC
%                InnEst.OptExitFlag
%                InnEst.OptOuPut

%  References
%   - "State and parameter estimation of stochastic physical systems 
%      from uncertain and indirect measurements",  
%      Jimenez J.C., Yoshimoto A., Miwakeichi F.
%      The European Physical Journal Plus, 
%      Volume 136, September 2021, 1-17.   
%   - "Bias reduction in the estimation of diffusion processes from discrete observations",
%      Jimenez J.C
%      IMA Journal of Mathematical Control and Information, 
%      Volumen 37, December 2020, Pages 1468–1505. 
%   - "Approximate linear minimum variance filters for continuous-discrete 
%      state space models: convergence and practical adaptive algorithms", 
%      Jimenez J.C.
%      IMA Journal of Mathematical Control and Information, 
%      Volume 36, Issue 2, June 2019, Pages 341–378. 
%
%  See also next reference for a review of inference methods for SDEs
%   - "Inference methods for discretely observed continuous-time stochastic volatility models: 
%     a commented overview”,
%     Jimenez J.C., Biscay R., Ozaki T., 
%     Asia-Pacific Financial Markets, 12 (2006) 109-141.
%

% The Likelihood function and the parameters to be estimated are normalized according the good optimization
% practice sugested at Chapter "Practicalities" in:
% Gill P.E., Murray W. and Wright M.H., Practical Optimization. Academic Press, 1981.

Xf00  = EstimOptions{3};
Theta0= EstimOptions{5};
LB = EstimOptions{10};
UB = EstimOptions{11};
OptimOptions = EstimOptions{12};
OptimMethod  = EstimOptions{17};
NumSim     = EstimOptions{20};

EstimOptions{13} = EstimChoise(ThetaIndex,Xf0Index);
EstimOptions{14} = Xf0Index;
EstimOptions{15} = ThetaIndex;

switch EstimOptions{13}
   case 0 
     n_par  = length(Xf0Index) + length(ThetaIndex);
     ParVar = [Xf00(Xf0Index)' Theta0(ThetaIndex)];
   case 1 
     n_par  = length(ThetaIndex);
     ParVar = Theta0(ThetaIndex);
   case 2 
     n_par  = length(Xf0Index);
     ParVar = Xf00(Xf0Index)';
end
NoZeroParVar = find(abs(ParVar)>eps);
ZeroParVar   = find(abs(ParVar)<eps);
  
if ~(isempty(LB) && isempty(UB))
  WrongLB=find(LB>ParVar);
  if ~isempty(WrongLB)
      disp('  ')
      disp(['A lower bound for the parameters is wrong. Please channge the bounds # : ',num2str(WrongLB)])
      ch='P';
      while ~(strcmp(ch,'Y') || strcmp(ch,'N'))
         ch=upper(input('The program will stop with error. Type "y" to continue ','s'));
      end
      return
  end 
  WrongUB=find(UB<ParVar);
  if ~isempty(WrongUB)
      disp('  ')
      disp(['A upper bound for the parameters is wrong. Please channge the bounds # : ',num2str(WrongUB)]) 
      ch='P';
      while ~(strcmp(ch,'Y') || strcmp(ch,'N'))
         ch=upper(input('The program will stop with error. Type "y" to continue ','s'));
      end
      return
  end 
  LB(NoZeroParVar) = LB(NoZeroParVar)./ParVar(NoZeroParVar);
  LB(ZeroParVar)   = LB(ZeroParVar) + 1;
  UB(NoZeroParVar) = UB(NoZeroParVar)./ParVar(NoZeroParVar);
  UB(ZeroParVar)   = UB(ZeroParVar) + 1;
  index=find(LB>UB);
  if ~isempty(index)
    temp=LB;
    LB(index)=UB(index);
    UB(index)=temp(index);
  end
end

figure
ParVar = ones(1,n_par);   
objfun = @(ParVar) LogLKFunction(ParVar,EstimOptions);
if isempty(EstimOptions{9})
  confun = [];
else     
  confun = @(ParVar) ConFunLike(ParVar,EstimOptions);
end    
warning('off','all');
switch OptimMethod
  case 'fmincon'
    [Par_est, ~, ExitFlag,OutPut] = fmincon(objfun,ParVar,[],[],[],[],LB,UB,confun,OptimOptions);
  case 'fminsearch'
    [Par_est, ~, ExitFlag,OutPut] = fminsearch(objfun,ParVar,OptimOptions);
  case 'patternsearch'
    [Par_est, ~, ExitFlag,OutPut] = patternsearch(objfun,ParVar,[],[],[],[],LB,UB,confun,OptimOptions);
  case 'UMDA+fmincon'
    if isempty(NumSim) || (length(NumSim)>1)
      [Par_est, ~] = GaussUMDA(@LogLKFunction,LB,UB,EstimOptions,OptimOptions);
    else
      [Par_est, ~] = GaussUMDA_parallel(@LogLKFunction,LB,UB,EstimOptions,OptimOptions);
    end
    [Par_est, ~, ExitFlag,OutPut] = fmincon(objfun,Par_est,[],[],[],[],LB,UB,confun,OptimOptions);
    OutPut.algorithm=['UMDA + ',OutPut.algorithm];
  otherwise
    error('Optimization Method is not well defined');    
end     
warning('on','all');
disp(['EXITFLAG: ',num2str(ExitFlag)])

t          = EstimOptions{1};
Z          = EstimOptions{2};
Pf0        = EstimOptions{4};
FileNames  = EstimOptions{6};
LKWindow   = EstimOptions{7};
mpts       = EstimOptions{16};
AbsTol     = EstimOptions{18};
RelTol     = EstimOptions{19};

Xf0_est=Xf00;
Theta_est=Theta0;  
switch EstimOptions{13}
  case 0
     NoZeroXf00   = intersect(find(abs(Xf00')>eps),Xf0Index);
     ZeroXf00     = intersect(find(abs(Xf00')<eps),Xf0Index);
     NoZeroTheta0 = intersect(find(abs(Theta0)>eps),ThetaIndex);
     ZeroTheta0   = intersect(find(abs(Theta0)<eps),ThetaIndex);
     AllPar=[Xf0_est' Theta_est];
     
     n_v_e=length(Xf00);
     AllPar([NoZeroXf00 n_v_e+NoZeroTheta0]) = [Xf00(NoZeroXf00)' Theta0(NoZeroTheta0)].* Par_est(NoZeroParVar);
     AllPar([ZeroXf00 ZeroTheta0])     = Par_est(ZeroParVar) - 1;  
     Xf0_est   = AllPar(1:n_v_e)';
     Theta_est = AllPar(1+n_v_e:end);
  case 1
     NoZeroTheta0 = intersect(find(abs(Theta0)>eps),ThetaIndex);
     ZeroTheta0   = intersect(find(abs(Theta0)<eps),ThetaIndex);
     Theta_est(NoZeroTheta0) = Theta0(NoZeroTheta0) .*  Par_est(NoZeroParVar);
     Theta_est(ZeroTheta0)   = Par_est(ZeroParVar) - 1;
  case 2
     NoZeroXf00 = intersect(find(abs(Xf00')>eps),Xf0Index);
     ZeroXf00   = intersect(find(abs(Xf00')<eps),Xf0Index);
     Xf0_est(NoZeroXf00) = Xf00(NoZeroXf00) .* Par_est(NoZeroParVar)';
     Xf0_est(ZeroXf00)   = Par_est(ZeroParVar)' - 1;
end

disp('Computing the Fitting Innovation')
if isempty(NumSim) || length(NumSim)>1
  LLF = LLfilter(t,Z,Xf0_est,Pf0,Theta_est,FileNames,LKWindow,mpts,AbsTol,RelTol,NumSim); 
elseif NumSim==0
  LLF = LLfilter(t,Z,Xf0_est,Pf0,Theta_est,FileNames,LKWindow,mpts);
else
  LLF = LLfilter(t,Z,Xf0_est,Pf0,Theta_est,FileNames,LKWindow,mpts,AbsTol,RelTol);
end

disp('Computing statistical tests for the Fitting Innovation')
alpha = 0.10; 
n_m = length(t);
NumPar = length(Xf0Index)+length(ThetaIndex);
for i=1:size(Z,1)
   obs= find(LLF.Inn(i,2:end)~=inf)+1;
   StdInn = LLF.Inn(i,obs)./sqrt(reshape(LLF.PInn(i,i,obs),1,length(obs)));
   [LLF.JBtest{i}.h,LLF.JBtest{i}.pValue,LLF.JBtest{i}.stat,LLF.JBtest{i}.cValue,Warning] = JBtest(StdInn',alpha);
   if Warning
     [LLF.JBtest{i}.h,LLF.JBtest{i}.pValue,LLF.JBtest{i}.stat,LLF.JBtest{i}.cValue] = JBtest(StdInn',alpha,0.0001);
   end
   [LLF.KStest{i}.h,LLF.KStest{i}.pValue,LLF.KStest{i}.stat,LLF.KStest{i}.cValue] = kstest(StdInn','alpha',alpha);
   MaxLag = maxLagEstimation(StdInn');
   LLF.LBQtest{i} = LBQtest(StdInn',MaxLag,NumPar,alpha);
   LLF.ARCHtest{i} = ARCHtest(StdInn','Lags',1:MaxLag,'Alpha',alpha);
end 

disp('Computing the Confidence Intervals')
[FMX,FMT] = FisherMatrix(t,Z,Xf0_est,Pf0,Theta_est,FileNames,LKWindow,mpts,AbsTol,RelTol,NumSim,Xf0Index,ThetaIndex,LLF.Inn,LLF.PInn);
Prob = 1 - alpha/2;
if ~isempty(FMX)
  FD = max([1 length(t)-length(Xf0Index)]);
  if FD==1, disp('Warning: number of parameters higher than number of observations'), end
  [U,D,~] = svd(FMX);
  D = 1./(diag(D)*FD);
  COV_Xf0 = U*diag(sqrt(D))*U';
  CI_Xf0  = tinv(Prob,FD-1)*diag(COV_Xf0);
else
  COV_Xf0 = [];
  CI_Xf0  = [];
end
if ~isempty(FMT)
  FD = max([1 length(t)-length(ThetaIndex)]);
  if FD==1, disp('Warning: number of parameters higher than number of observations'), end
  [U,D,~] = svd(FMT);
  D = 1./(diag(D)*FD);
  COV_Theta = U*diag(sqrt(D))*U';
  CI_Theta  = tinv(Prob,FD-1)*diag(COV_Theta)';
else
  COV_Theta = [];
  CI_Theta  = [];
end
AIC = 2*NumPar + LLF.LogLK;
BIC = NumPar.*log(n_m) + LLF.LogLK;

if nargout==1
    varargout{1}.Xf0      = Xf0_est;
    varargout{1}.Theta    = Theta_est;
    varargout{1}.LLF      = LLF;
    varargout{1}.Xf0CI    = CI_Xf0;
    varargout{1}.ThetaCI  = CI_Theta;
    varargout{1}.Xf0Cov   = COV_Xf0;
    varargout{1}.ThetaCov = COV_Theta;
    varargout{1}.AIC      = AIC;
    varargout{1}.BIC      = BIC;
    varargout{1}.OptExitFlag = ExitFlag;
    varargout{1}.OptOutPut   = OutPut;
   
else
    varargout{1} = Xf0_est;
    varargout{2} = Theta_est;
    varargout{3} = LLF;
    varargout{4} = CI_Xf0;
    varargout{5} = CI_Theta;
    varargout{6} = COV_Xf0;
    varargout{7} = COV_Theta;
    varargout{8} = AIC;
    varargout{9} = BIC;
end



%-------------------------------------------

function option = EstimChoise(ThetaIndex,Xf0Index)

if (~isempty(ThetaIndex)) && (isempty(Xf0Index))  
   option=1;
else  
   if (isempty(ThetaIndex)) && (~isempty(Xf0Index)) 
      option=2;
   else
      option=0;
   end
end  




%----------------------------
% Computation of the Fisher Information Matrix for the estimated parameters 

function V = auxiliar(t,ParVar,Xf0_est,Xf0Index,Theta_est,ThetaIndex,FilterType)

option = ~isempty(Xf0Index) - ~isempty(ThetaIndex); 
switch option
   case 0 
     Xf0_est(Xf0Index) = ParVar(1:length(Xf0Index));
     Theta_est(ThetaIndex) = ParVar(length(Xf0Index)+1:end);
     LLF = FilterType(Xf0_est,Theta_est');
   case -1 
     Theta_est(ThetaIndex) = ParVar;
     LLF = FilterType(Theta_est');
   case 1 
     Xf0_est(Xf0Index) = ParVar;
     LLF = FilterType(Xf0_est);
end

[NV,NO]=size(LLF.Inn);
V = [reshape(LLF.Inn,NV*NO,1); reshape(LLF.PInn,NV*NV*NO,1)];

                        
function [FMX,FMT] = FisherMatrix(t,Z,Xf0_est,Pf0,Theta_est,FileNames,LKWindow,mpts,AbsTol,RelTol,NumSim,Xf0Index,ThetaIndex,FInn,PFInn)

option = ~isempty(Xf0Index) - ~isempty(ThetaIndex); 

switch option
    case 0
        ParVar = [Xf0_est(Xf0Index); Theta_est(ThetaIndex)'];
        if isempty(NumSim) || length(NumSim)>1
            FilterType = @(Xf0_est,Theta_est) LLfilter(t,Z,Xf0_est,Pf0,Theta_est,FileNames,LKWindow,mpts,AbsTol,RelTol,NumSim); 
          elseif NumSim==0
             FilterType = @(Xf0_est,Theta_est) LLfilter(t,Z,Xf0_est,Pf0,Theta_est,FileNames,LKWindow,mpts);
          else
             FilterType = @(Xf0_est,Theta_est) LLfilter(t,Z,Xf0_est,Pf0,Theta_est,FileNames,LKWindow,mpts,AbsTol,RelTol);
        end
    case -1
        ParVar = Theta_est(ThetaIndex)';
        if isempty(NumSim) || length(NumSim)>1
            FilterType = @(Theta_est) LLfilter(t,Z,Xf0_est,Pf0,Theta_est,FileNames,LKWindow,mpts,AbsTol,RelTol,NumSim); 
          elseif NumSim==0
             FilterType = @(Theta_est) LLfilter(t,Z,Xf0_est,Pf0,Theta_est,FileNames,LKWindow,mpts);
          else
             FilterType = @(Theta_est) LLfilter(t,Z,Xf0_est,Pf0,Theta_est,FileNames,LKWindow,mpts,AbsTol,RelTol);
        end
    case 1
        ParVar = Xf0_est(Xf0Index);
        if isempty(NumSim) || length(NumSim)>1
            FilterType = @(Xf0_est) LLfilter(t,Z,Xf0_est,Pf0,Theta_est,FileNames,LKWindow,mpts,AbsTol,RelTol,NumSim); 
          elseif NumSim==0
             FilterType = @(Xf0_est) LLfilter(t,Z,Xf0_est,Pf0,Theta_est,FileNames,LKWindow,mpts);
          else
             FilterType = @(Xf0_est) LLfilter(t,Z,Xf0_est,Pf0,Theta_est,FileNames,LKWindow,mpts,AbsTol,RelTol);
        end
end

[NV,NO]=size(FInn);
V0 = [reshape(FInn,NV*NO,1); reshape(PFInn,NV*NV*NO,1)];
THRESH=10^-9*ones(length(ParVar),1);
funcion=@(tt,ParVar) auxiliar(tt,ParVar,Xf0_est,Xf0Index,Theta_est,ThetaIndex,FilterType);
[DFDY,~] = numjac(funcion,[],ParVar,V0,THRESH,[],0);

NXf0=length(Xf0Index);
NTheta=length(ThetaIndex);
switch option
     case 0       
        JmX=DFDY(1:NV*NO,1:NXf0);
        JmX=reshape(JmX,NV,NO,NXf0);
        JPX=DFDY(NV*NO+1:end,1:NXf0);
        JPX=reshape(JPX,NV*NV,NO,NXf0);
        
        JmT=DFDY(1:NV*NO,NXf0+1:end);
        JmT=reshape(JmT,NV,NO,NTheta);
        JPT=DFDY(NV*NO+1:end,NXf0+1:end);
        JPT=reshape(JPT,NV*NV,NO,NTheta);
     case -1             
        JmT=DFDY(1:NV*NO,:);
        JmT=reshape(JmT,NV,NO,NTheta);
        JPT=DFDY(NV*NO+1:end,:);
        JPT=reshape(JPT,NV*NV,NO,NTheta);
     case 1
        JmX=DFDY(1:NV*NO,:);
        JmX=reshape(JmX,NV,NO,NXf0);
        JPX=DFDY(NV*NO+1:end,:);
        JPX=reshape(JPX,NV*NV,NO,NXf0);
end  

if option == 0 || option == -1
  FMT=zeros(NTheta);
  for k=2:NO
     obs = find(FInn(:,k)~=inf);
%     NV = length(obs);
     Inv = inv(PFInn(obs,obs,k));
     for i=1:NTheta
        VI = JmT(obs,k,i)' * Inv;
        MI = reshape(JPT(:,k,i),NV,NV);
        MI = Inv * MI(obs,obs); 
        for j=i:NTheta
           MJ = reshape(JPT(:,k,j),NV,NV);
           MJ = Inv * MJ(obs,obs); 
           FMT(i,j) = FMT(i,j) + VI*JmT(obs,k,j) + 0.5*trace(MI*MJ);
         end
      end
   end
  FMT = FMT + FMT' -diag(diag(FMT));
else
  FMT=[];
end

if option == 0 || option == 1
  FMX=zeros(NXf0);
  for k=2:NO
     obs = find(FInn(:,k)~=inf);
     NV = length(obs);
     Inv = inv(PFInn(obs,obs,k));
     for i=1:NXf0
        VI = JmX(obs,k,i)' * Inv;
        MI = reshape(JPX(:,k,i),NV,NV);
        MI = Inv * MI(obs,obs); 
       for j=i:NXf0
           MJ = reshape(JPX(:,k,j),NV,NV);
           MJ = Inv * MJ(obs,obs); 
           FMX(i,j) = FMX(i,j) + VI*JmX(obs,k,j) + 0.5*trace(MI*MJ);
        end
     end
  end
  FMX = FMX + FMX' -diag(diag(FMX));
else
  FMX=[];
end



function maxLag = maxLagEstimation(res)

%   [1] Box, G.E.P., G.M. Jenkins, and G.C. Reinsel. Time Series Analysis:
%       Forecasting and Control. 3rd ed. Upper Saddle River, NJ:
%       Prentice-Hall, 1994.

h=0;
maxLag=max([1 min([20 length(res)-3])]);   % 20 is recomended in [1]. Following archtest.m, maxLag must be less than length(res)-2
specU = garch(0,maxLag);
[~,~,uLL] = estimate(specU,res,'Display','off');
while h==0 && maxLag>1
  maxLag=maxLag-1;
  specR = garch(0,maxLag);
  [~,~,rLL] = estimate(specR,res,'Display','off');
  if rLL < uLL
    h = lratiotest(uLL,rLL,1);
  end
  uLL = rLL;
end



  