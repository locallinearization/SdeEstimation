function PlotModelEstimates(t,EF,LLF,MomentType,var)
%
% Plots exact and approximate conditional moments
%
% Input variables
%           t : time
%          EF : exact filter
%         LLF : approximate filter
%  MomentType : string with two options 'First' or 'Second' conditional moment
%         var : n_plots x 1 vector with the numbers of the state variables to plot
%                 or
%         var : n_plots x 2 matrix with the entries of the covariance variables to plot
%  
% Call options
%   PlotModelEstimates(t,EF,LLF,'First') 
%          plot the exact and LL filter, and the exact and LL prediction for all the state variables
%  
%   PlotModelEstimates(t,EF,LLF,'First',var)
%          plot the exact and LL filter, and the exact and LL prediction for the state variables
%          specified in vector var
%
%   PlotModelEstimates(t,EF,LLF,'Second')
%          plot the diagonal of the exact and LL filter variance, and the diagonal of the exact and LL prediction variance 
%
%   PlotModelEstimates(t,EF,LLF,'Second',var)
%          plot the exact and LL filter covariance entries, and the exact and LL prediction covariance entries 
%          specified in the matrix var

switch MomentType
    case 'First'
       switch nargin
           case 4
             NumSubplots=size(LLF.Xf,1);
             var=1:NumSubplots;
           case 5
             NumSubplots=length(var);
       end
    case 'Second'
       switch nargin
           case 4
             NumSubplots=size(LLF.Xf,1);
             var=1:NumSubplots;
             Row=var;
             Col=var;
           case 5
             NumSubplots=size(var,1);
             Row=var(:,1);
             Col=var(:,2);
             var=1:NumSubplots;
       end
end

if strcmp(MomentType,'First')
    n=0;
    figure
    for i=var
        n = n+1;
        Str = [num2str(NumSubplots) '1' num2str(n)];
        subplot(Str)    
        plot(t,EF.Xf(i,:),t,LLF.Xf(i,:));
        ylabel(['x',num2str(i)])
        if isfield(LLF,'JBtest')
          title('exact (blue) vs. LL fitting-filter')
        else
          title('exact (blue) vs. LL filter')
        end
    end
    xlabel('t')

    n=0;
    figure
    for i=var
        n = n+1;
        Str = [num2str(NumSubplots) '1' num2str(n)];
        subplot(Str)    
        plot(t,EF.Xp(i,:),t,LLF.Xp(i,:));
        ylabel(['x',num2str(i)])
        if isfield(LLF,'JBtest')
          title('exact (blue) vs. LL fitting-prediction')
        else
          title('exact (blue) vs. LL prediction')
        end

    end
    xlabel('t')
else
    figure
    n_m = length(t);
    for i=var
        Str = [num2str(NumSubplots) '1' num2str(i)];
        subplot(Str)    
        v1=reshape(EF.Pf(Row(i),Col(i),:),1,n_m);
        v2=reshape(LLF.Pf(Row(i),Col(i),:),1,n_m);
        plot(t,v1,t,v2);
        ylabel(['P',num2str(Row(i)),num2str(Col(i))])
        if isfield(LLF,'JBtest')
          title('exact (blue) vs. LL fitting-filter variance')
        else
          title('exact (blue) vs. LL filter variance')
        end
    end
    xlabel('t')

    figure
    for i=var
        Str = [num2str(NumSubplots) '1' num2str(i)];
        subplot(Str)    
        v1=reshape(EF.Pp(Row(i),Col(i),:),1,n_m);
        v2=reshape(LLF.Pp(Row(i),Col(i),:),1,n_m);
        plot(t,v1,t,v2);
        ylabel(['P',num2str(Row(i)),num2str(Col(i))])
        if isfield(LLF,'JBtest')
          title('exact (blue) vs. LL fitting-prediction variance')
        else
          title('exact (blue) vs. LL prediction variance')
        end
    end
    xlabel('t')
end
