function PlotFilterEstimates(t,Z,LLF,IND_OBS_X)
%
% Plots the observations Z and their filter estimates LLF.
%         t : time
%         Z : observed data
%       LLF : filtering estimates
% IND_OBS_X : Vector of boolean variables defining whether each vector
%             component of the state variable have been observed or not

if isfield(LLF,'AlertPts') && LLF.AlertPts
   disp('  ')
   disp('warning: accepted + rejected points reachs the allowed maximum')
end 

if isfield(LLF,'AlertSim') && LLF.AlertSim
   disp('  ')
   disp('warning: number of simulations reachs the allowed maximum')
end 
disp('  ')

n_obs = sum(IND_OBS_X);
NZ=find(IND_OBS_X>0);
NumSubplots=size(LLF.pts,1)+2;
if n_obs ~= size(Z,1)
    error('Length of vector "IND_OBS_X" wrong')
end

% plot of the state variable and their filter estimates
for i=1:n_obs
    figure
    Str = [num2str(NumSubplots) '11'];
    subplot(Str)    
    obs= find(Z(i,:)~=inf);
    plot(t(obs),Z(i,obs),t(obs),LLF.Xf(NZ(i),obs));
    ejes=axis;
    axis([t(1) t(end) ejes(3:4)])
    ylabel(['z',num2str(i)])
    if isfield(LLF,'JBtest')
      title('data (blue) vs. LL fitting-filter')
    else
      title('data (blue) vs. LL filter')
    end
    
    Str = [num2str(NumSubplots) '12'];
    subplot(Str)
    obs= find(LLF.Inn(i,:)~=inf);
    plot(t(obs),LLF.Inn(i,obs));
    ejes=axis;
    axis([t(1) t(end) ejes(3:4)])
    subplot(Str) 
    if isfield(LLF,'JBtest')
       title('fitting-innovation')
    else
       title('innovation')
    end
    
    Str = [num2str(NumSubplots) '13'];
    subplot(Str)
    plot(t,LLF.pts(1,:))
    ejes=axis;
    axis([t(1) t(end) ejes(3:4)])
    Min=num2str(min(LLF.pts(1,2:end)));
    Max=num2str(max(LLF.pts(1,2:end)));
    Sum=num2str(sum(LLF.pts(1,:)));
    title(['apts,     Min= ',Min,'   Max= ',Max,'  Total= ',Sum])
    
    if size(LLF.pts,1)>1 
        Str = [num2str(NumSubplots) '14'];
        subplot(Str)
        plot(t,LLF.pts(2,:))
        ejes=axis;
        axis([t(1) t(end) ejes(3:4)])
        MinF=num2str(min(LLF.pts(2,2:end)));
        MaxF=num2str(max(LLF.pts(2,2:end)));
        SumF=num2str(sum(LLF.pts(2,:)));
        title(['fpts,     Min= ',MinF,'   Max= ',MaxF,'  Total= ',SumF])

        if size(LLF.pts,1)>2
            Str = [num2str(NumSubplots) '15'];
            subplot(Str)
            plot(t,LLF.pts(3,:))
            ejes=axis;
            axis([t(1) t(end) ejes(3:4)])
            MinS=num2str(min(LLF.pts(3,2:end)));
            MaxS=num2str(max(LLF.pts(3,2:end)));
            SumS=num2str(sum(LLF.pts(3,:)));
            title(['NSim,     Min= ',MinS,'   Max= ',MaxS,'  Total= ',SumS])
        end
    end
    xlabel('t')
end

