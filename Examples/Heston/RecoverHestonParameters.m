% script for getting back the estimation of the original Heston parameter values

B1=InnEst.Theta(2);
L=InnEst.Theta(3);         % lamda
M=B1./L;                   % mu
dB1=InnEst.ThetaCI(2);
dL=InnEst.ThetaCI(3);
dM=(dB1-M.*dL)./L;

B2=InnEst.Theta(4);
B3=InnEst.Theta(5);

G=sqrt(B2.*B2 + B3.*B3);   % gamma
R=B2./G;                   % ro

dB2=InnEst.ThetaCI(4);
dB3=InnEst.ThetaCI(5);
c = sqrt(1-R*R);
Delta = [R G; c -G*R/c ]\[dB2; dB3];

HestonInnEst.Theta=InnEst.Theta;
HestonInnEst.Theta([2 4 5])=[M G R];
HestonInnEst.ThetaCI=InnEst.ThetaCI;
HestonInnEst.ThetaCI([2 4 5])=abs([dM Delta']);

HestonInnEst.LLF.LogLK=InnEst.LLF.LogLK;
HestonInnEst.AIC=InnEst.AIC;
HestonInnEst.BIC=InnEst.BIC;
HestonTheta=[alpha mu lamda gamma ro];
EstimOptions{17}='UMDA+fmincon';

disp([' Original Heston parameter values'])
DisplayParameterEstimates(HestonInnEst,X0,HestonTheta,Xf0Index,ThetaIndex,EstimOptions)

