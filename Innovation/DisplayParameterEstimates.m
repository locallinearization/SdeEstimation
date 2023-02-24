function DisplayParameterEstimates(InnEst,X0,Theta,Xf0Index,ThetaIndex,EstimOptions)
%
% Summarizes the parameter estimation 
%
% Input
%    InnEst:       output of the function "InnovEstimator"
%    X0:           initial state variable
%    Theta:        values of the state space parameters
%    Xf0Index:     indexes of the estimated initial state variables
%    ThetaIndex:   indexes of the estimated parameters
%    EstimOptions: output of the function "EstimSet" used to obtain InnEst 

Xf0=EstimOptions{3};
Theta0=EstimOptions{5};
OptimMethod=EstimOptions{17};

disp(['   '])
disp([' Estimation Summary  '])
disp(['   '])
disp([' Log LK = ' num2str(-0.5*InnEst.LLF.LogLK)])
disp(['    AIC = ' num2str(InnEst.AIC)])
disp(['    BIC = ' num2str(InnEst.BIC)])
disp(['   '])
if ~isempty(ThetaIndex)
  ST.ParNum     = ThetaIndex';
  ST.Theta      = Theta(ThetaIndex)';
  if ~strcmp(OptimMethod,'UMDA+fmincon')
    ST.Theta0     = Theta0(ThetaIndex)';
  end
  ST.Theta_Est  = InnEst.Theta(ThetaIndex)';
  ST.Conf_Inter = InnEst.ThetaCI';
  TableT = struct2table(ST);
  disp(TableT);
end
if ~isempty(Xf0Index)
  SX.StateNum   = Xf0Index';
  SX.X          = X0(Xf0Index);
  if ~strcmp(OptimMethod,'UMDA+fmincon')
    SX.X0         = Xf0(Xf0Index);
  end
  SX.X_Est      = InnEst.Xf0(Xf0Index);
  SX.Conf_Inter = InnEst.Xf0CI;
  TableX = struct2table(SX);
  disp(TableX);
end
