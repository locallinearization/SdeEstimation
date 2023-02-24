function LogLK = LogLKFunction(ThetaVar,EstimOptions)

% The Likelihood function and the parameters to be estimated are normalized according the good optimization
% practice sugested at Chapter "Practicalities" in:
% Gill P.E., Murray W. and Wright M.H., Practical Optimization. Academic Press, 1981.

global gct    % variable to be used in confun.m

t          = EstimOptions{1};
Z          = EstimOptions{2};
Xf00       = EstimOptions{3};
Pf0        = EstimOptions{4};
Theta0     = EstimOptions{5};
FileNames  = EstimOptions{6};
LKWindow   = EstimOptions{7};
LogLK0     = EstimOptions{8};
nu_stast   = EstimOptions{9};
Xf0Index   = EstimOptions{14};
ThetaIndex = EstimOptions{15};
mpts       = EstimOptions{16};
AbsTol     = EstimOptions{18};
RelTol     = EstimOptions{19};
NumSim     = EstimOptions{20};

switch EstimOptions{13}
   case 0
      Xf  = Xf00;
      Th  = Theta0;
      Par = [Xf00(Xf0Index)' Theta0(ThetaIndex)];  
      NoZeroPar = find(abs(Par)>eps);
        ZeroPar = find(abs(Par)<eps);
        
      Par(NoZeroPar)= ThetaVar(NoZeroPar) .* Par(NoZeroPar);
      Par(ZeroPar)  = ThetaVar(ZeroPar) - 1;
      Xf(Xf0Index)  = Par(1:length(Xf0Index));
      Th(ThetaIndex)= Par(1+length(Xf0Index):end);
      if isempty(NumSim) || length(NumSim)>1
         LLF=LLfilter(t,Z,Xf,Pf0,Th,FileNames,LKWindow,mpts,AbsTol,RelTol,NumSim); 
      elseif NumSim==0
          LLF=LLfilter(t,Z,Xf,Pf0,Th,FileNames,LKWindow,mpts);
      else
          LLF=LLfilter(t,Z,Xf,Pf0,Th,FileNames,LKWindow,mpts,AbsTol,RelTol);
      end
   case 1
      Th  = Theta0;
      Par = Theta0(ThetaIndex);         % cada elemento de Par tiene su correspondiente elemento en ThetaVar
      NoZeroPar = find(abs(Par)>eps);
        ZeroPar = find(abs(Par)<eps);

      Par(NoZeroPar)= ThetaVar(NoZeroPar) .* Par(NoZeroPar);
      Par(ZeroPar)  = ThetaVar(ZeroPar) - 1;
      Th(ThetaIndex)= Par;
      if isempty(NumSim) || length(NumSim)>1
         LLF=LLfilter(t,Z,Xf00,Pf0,Th,FileNames,LKWindow,mpts,AbsTol,RelTol,NumSim); 
      elseif NumSim==0
          LLF=LLfilter(t,Z,Xf00,Pf0,Th,FileNames,LKWindow,mpts);
      else
          LLF=LLfilter(t,Z,Xf00,Pf0,Th,FileNames,LKWindow,mpts,AbsTol,RelTol);
      end
   case 2
      Xf  = Xf00;
      Par = Xf00(Xf0Index);            % cada elemento de Par tiene su correspondiente elemento en ThetaVar
      NoZeroPar = find(abs(Par)>eps);
        ZeroPar = find(abs(Par)<eps);
      
      Par(NoZeroPar)= ThetaVar(NoZeroPar)' .* Par(NoZeroPar);
      Par(ZeroPar)  = ThetaVar(ZeroPar) - 1;
      Xf(Xf0Index)= Par;
      if isempty(NumSim) || length(NumSim)>1
         LLF=LLfilter(t,Z,Xf,Pf0,Theta0,FileNames,LKWindow,mpts,AbsTol,RelTol,NumSim); 
      elseif NumSim==0
          LLF=LLfilter(t,Z,Xf,Pf0,Theta0,FileNames,LKWindow,mpts);
      else
          LLF=LLfilter(t,Z,Xf,Pf0,Theta0,FileNames,LKWindow,mpts,AbsTol,RelTol);
      end
    otherwise
end
   
PlotStandInn(t,LLF.Inn,LLF.PInn,LLF.pts,LKWindow)

if ~isempty(nu_stast)
  % see in matlab help "Passing Extra Parameters" 
  index=find((t>=LKWindow(1)) & (t<=LKWindow(2)));
  firstpoint=index(1);                 % starting point for the computation of the innovations statistics
  lastpoint=index(end);                % last point for the computation of the innovations statistics
  nut=LLF.Inn(:,firstpoint:lastpoint)';
  gct=[abs(mean(nut)) diag(cov(nut))'];
  gct=gct-nu_stast;
end

if (exist('LLF.AlertPts','var') && LLF.AlertPts) || (exist('LLF.AlertSim','var') && LLF.AlertSim)
    LogLK = NaN;
else
    LogLK = LLF.LogLK./abs(LogLK0);
end



function PlotStandInn(t,Inn,PInn,pts,LKWindow)

n_obs = size(Inn,1);
NumSubplots=n_obs+size(pts,1);
for i=1:n_obs   
   Str = [num2str(NumSubplots) '1' num2str(i)];
   subplot(Str)
   obs= find(Inn(i,:)~=inf);
   plot(t(obs),Inn(i,obs)./sqrt(reshape(PInn(i,i,obs),1,length(t(obs)))))
   NoLK = find((t(obs)<LKWindow(1)) | (t(obs)>LKWindow(2)) );
   if ~isempty(NoLK)
       hold on 
       plot(t(NoLK),Inn(i,NoLK)./sqrt(reshape(PInn(i,i,NoLK),1,length(t(NoLK)))),'.r')
       hold off
   end  
   ejes=axis;
   axis([t(1) t(end) ejes(3:4)])
   title('Standardized fitting-innovation inside the optimization algorithm')
end

Str = [num2str(NumSubplots) '1' num2str(n_obs+1)];
subplot(Str)
plot(t,pts(1,:))
ejes=axis;
axis([t(1) t(end) ejes(3:4)])
Min=num2str(min(pts(1,2:end)));
Max=num2str(max(pts(1,2:end)));
Sum=num2str(sum(pts(1,:)));
title(['apts,     Min= ',Min,'   Max= ',Max,'  Total= ',Sum])
    
if size(pts,1)>1 
   Str = [num2str(NumSubplots) '1' num2str(n_obs+2)];
   subplot(Str)
   plot(t,pts(2,:))
   ejes=axis;
   axis([t(1) t(end) ejes(3:4)])
   MinF=num2str(min(pts(2,2:end)));
   MaxF=num2str(max(pts(2,2:end)));
   SumF=num2str(sum(pts(2,:)));
   title(['fpts,     Min= ',MinF,'   Max= ',MaxF,'  Total= ',SumF])

   if size(pts,1)>2
      Str = [num2str(NumSubplots) '1' num2str(n_obs+3)];
      subplot(Str)
      plot(t,pts(3,:))
      ejes=axis;
      axis([t(1) t(end) ejes(3:4)])
      MinS=num2str(min(pts(3,2:end)));
      MaxS=num2str(max(pts(3,2:end)));
      SumS=num2str(sum(pts(3,:)));
      title(['NSim,     Min= ',MinS,'   Max= ',MaxS,'  Total= ',SumS])
    end
end
xlabel('t')
drawnow



