function [Zt, Vt, nu, Sv] = Filter_Add(tt,Ztdt,Pt,theta,X,FileName)

% LL filter estimates for Observations with Additive Noise

[he,Jhe,See] = FileName(tt,Ztdt,theta);

obs= find(X(:)~=inf);
he = he(obs);
Jhe = Jhe(obs,:);
See = See(obs,obs);

nu = X(obs) - he;           %  innovation
Sv   = Jhe*Pt*Jhe' + See;   %  innovation variance             

Kt = Pt*Jhe' / Sv;                      %(4)    LL filter gain
Zt = Ztdt + Kt * nu;                    %(5)    LL filter
Vt = Pt - Kt * Jhe*Pt;                  %(6)    LL filter variance

