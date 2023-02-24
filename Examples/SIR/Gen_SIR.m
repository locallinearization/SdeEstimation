% Script for model simulation
close all;

% Model definition
Const_SIR;

% observed noise seed
ONseed=100; 

% integration by RK method
tspan=t0:deltat:tTotal;
ode=@SE_SIR;
[t,proceso]=ode45(ode,tspan,X0,[],Theta);

% sampling the observations for filtering
t=t';
proceso=proceso';

% generation of noisy observations
datos=Observations(t,proceso,Theta,FileNames,ObsX,ONseed);   

% selecting the missing data
n_m=length(t);
v=rand(1,n_m);
index=find(v>0.5);
datos(1,index(2:end))=inf;    % setting the missing data for the observations of the state variable 'I'
index=find(v<0.1);
datos(2,index(2:end))=inf;    % setting the missing data for the observations of the state variable 'R'

% plot state variables and observations
PlotModelRealization(t,datos,proceso,ObsX);

% save generated data
FilePath=fileparts(which(ModEq));
Name=[FilePath,filesep,'data' ModEq(3:end)];
save(Name,'t','datos','proceso') 
