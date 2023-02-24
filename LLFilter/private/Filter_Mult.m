function [Zt, Vt, nu, Sv] = Filter_Mult(tt,Ztdt,Pt,theta,X,FileName)

% LL filter estimates for Observations with Multiplicative Noise

[he,Jhe,See,pe,Jpe,Lee] = FileName(tt,Ztdt,theta);

obs= find(X(:)~=inf);
he = he(obs);
Jhe = Jhe(obs,:);
See = See(obs,obs);

nu = X(obs) - he;           %  innovation

n=size(pe,2);
p_i=See;
for j=1:n
   pe{j} = pe{j}(obs);
   Jpe{j} = Jpe{j}(obs,:);
   Lee{j} = Lee{j}(obs,obs);

   th{j} = pe{j}-Jpe{j}*Ztdt;
   p_i   = p_i + Lee(j)*th{j}*th{j}';
   th{j} = Lee(j)*th{j};
end
Sv   = Jhe*Pt*Jhe' + p_i;                  % es parte de (4)
for j=1:n
   temp=Jpe{j}*Ztdt*th{j}';
   Sv = Sv + Lee(j)*Jpe{j}*(Pt+Ztdt*Ztdt')*Jpe{j}' + temp + temp';     % innovation variance
end  
Kt = Pt*Jhe' / Sv;                      %(4)    LL filter gain
Zt = Ztdt + Kt * nu;                    %(5)    LL filter
Vt = Pt - Kt * Jhe*Pt;                  %(6)    LL filter variance

