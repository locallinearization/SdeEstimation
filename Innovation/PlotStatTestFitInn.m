function PlotStatTestFitInn(LLF)
%
% plot the statistics of the standarized fitting-innovation:
%
%      Jarque-Bera test of composite Gaussianity
%      Kolmogorov-Smirnov test of normality N(0,1)
%      Ljung-Box Q-test for the autocorrelation 
%      Ljung-Box Q-test of heteroscedasticity

%  Reference
%   - "State and parameter estimation of stochastic physical systems 
%      from uncertain and indirect measurements",  
%      Jimenez J.C., Yoshimoto A., Miwakeichi F.
%      The European Physical Journal Plus, 
%      Volume 136, September 2021, 1-17.   


[n_obs,~] = size(LLF.Inn);
for i=1:n_obs
    figure
    subplot(3,1,1)
    obs= find(LLF.Inn(i,2:end)~=inf)+1;
    StdInn = LLF.Inn(i,obs)./sqrt(reshape(LLF.PInn(i,i,obs),1,length(obs)));
    hist(StdInn)
    if LLF.JBtest{i}.h, resultJB='rejected'; else resultJB='accepted'; end
    if LLF.KStest{i}.h, resultKS='rejected'; else resultKS='accepted'; end
    if n_obs>1
      title(['Histogram of the standarized fitted-innovation ',num2str(i)])
    else
      title('Histogram of the standarized fitted-innovation ')
    end
    stringJBtest=['Gaussianity test = ',resultJB,', cValue = ',num2str(LLF.JBtest{i}.cValue),', stat = ',num2str(LLF.JBtest{i}.stat)];
    stringKStest=['Kolmogorov-Smirnov test = ',resultKS,', cValue = ',num2str(LLF.KStest{i}.cValue),', stat = ',num2str(LLF.KStest{i}.stat)];
    xlabel({stringJBtest,stringKStest})

    rejected=LLF.LBQtest{i}.h;
    accepted=~LLF.LBQtest{i}.h;
    xpoints=1:length(rejected);
    subplot(3,1,2)
    plot(LLF.LBQtest{i}.cValue,'s-k')
    hold on
    caso = (sum(rejected)>0) - (sum(accepted)>0); 
    switch caso
        case -1
           plot(xpoints(accepted),LLF.LBQtest{i}.stat(accepted),'*b')
           legend('cValues','Accepted','Location','Best')
        case  0
           plot(xpoints(accepted),LLF.LBQtest{i}.stat(accepted),'*b')
           plot(xpoints(rejected),LLF.LBQtest{i}.stat(rejected),'*r')
           legend('cValues','Accepted','Rejected','Location','Best')
       case  1
           plot(xpoints(rejected),LLF.LBQtest{i}.stat(rejected),'*r')
           legend('cValues','Rejected','Location','Best')
    end
    xlabel('lag')
    ylabel('stat')
    title('LBQ test')

    rejected=LLF.ARCHtest{i}.h;
    accepted=~LLF.ARCHtest{i}.h;
    subplot(3,1,3)
    plot(LLF.ARCHtest{i}.cValue,'s-k')
    hold on
    caso = (sum(rejected)>0) - (sum(accepted)>0); 
    switch caso
        case -1
           plot(xpoints(accepted),LLF.ARCHtest{i}.stat(accepted),'*b')
           legend('cValues','Accepted','Location','Best')
        case  0
           plot(xpoints(accepted),LLF.ARCHtest{i}.stat(accepted),'*b')
           plot(xpoints(rejected),LLF.ARCHtest{i}.stat(rejected),'*r')
           legend('cValues','Accepted','Rejected','Location','Best')
       case  1
           plot(xpoints(rejected),LLF.ARCHtest{i}.stat(rejected),'*r')
           legend('cValues','Rejected','Location','Best')
    end
    xlabel('lag')
    ylabel('stat')
    title('ARCH test')
end