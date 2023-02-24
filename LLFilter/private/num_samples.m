function [NS_out, NoEnd, Info]=num_samples(M1,M2,M4,NS_in,SRTol,SATol,NumWorkers,MaxNumSim)
% Estimation of a new number of stochatic simulations
% 
% Input
%   M1,M2,M4: first, second and fourth moments
%   NS_in: current number of simulation
%   SRTol,SATol: relative and absolute error for the estimation of the number of simulations
%   NumWorkers: number of workers
%   MaxNumSim: maximum number of allowed simulations
%
% Output
%   NS_out: new number of stochatic simulations 
%   NoEnd : binary variable (0: end of the stochastic simulations)
%   Info : Information about the statistical errors

% Comm Pure & Appl. Math., LIV (2001) 1169-1214, but estimating the final number 
% of simulations by relative error criteria (instead of the Szepessy's absolute error criteria 
% used in IMA J.Control (2019)).

c0=1.65;

dMed=abs(M2-M1.*M1);
dVar=abs(M4-M2.*M2);
ErrorMed=c0*sqrt(dMed./(NS_in*(SATol + SRTol*abs(M1.*M1))));
ErrorVar=c0*sqrt(dVar./(NS_in*(SATol + SRTol*abs(M2.*M2))));
NoEndErrorMed=max( ErrorMed )>=1;
NoEndErrorVar=max( ErrorVar )>=1;
NoEnd = NoEndErrorMed || NoEndErrorVar;
NoEnd = NoEnd & NS_in~=MaxNumSim;

BS='           ';
if NoEnd
  cteMed=(c0/(0.95*(SATol + SRTol*abs(M1.*M1)))).^2; 
  cteVar=(c0/(0.95*(SATol + SRTol*abs(M2.*M2)))).^2;
  NS_Med=min([max(round(cteMed*dMed)); 2*NS_in]);
  NS_Var=min([max(round(cteVar*dVar)); 2*NS_in]);
  NS_out=max([NS_Med; NS_Var; NS_in+NumWorkers]);
  NS_out=min([NS_out MaxNumSim]);
  if ~NoEndErrorMed, NS_Med=0; end
  if ~NoEndErrorVar, NS_Var=0; end
  Info=[num2str(max(ErrorMed)),BS,num2str(max(ErrorVar)),BS,num2str(NS_Med),BS,num2str(NS_Var),BS,num2str(NS_out-NS_in)];
else
  NS_out=NS_in;
  Info=[num2str(max(ErrorMed)),BS,num2str(max(ErrorVar)),BS,num2str(0),BS,num2str(0),BS,num2str(NS_out)];
end





