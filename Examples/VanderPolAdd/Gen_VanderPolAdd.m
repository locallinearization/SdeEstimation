% Script for model simulation
close all

% Model definition
Const_VanderPolAdd;

% system noise seed
SNseed=1;  
% observed noise seed
ONseed=31; 

% number of points for filtering
n_m=round((tTotal - t0) ./ deltat) + 1;

% number of points for simulations
n_m_max=round((tTotal - t0) ./ deltat_min) + 1;
pp=(n_m_max-1)/(n_m-1);

% integration by Euler method
[t,proceso]=Euler(t0,X0,n_m_max,deltat_min,FileNames,Theta,SNseed); 
% sampling the observations for filtering
t=t(1:pp:n_m_max);                    
proceso=proceso(:,1:pp:n_m_max);

% generation of noisy observations
datos=Observations(t,proceso,Theta,FileNames,ObsX,ONseed);   

% plot state variables and observations
PlotModelRealization(t,datos,proceso,ObsX);

% save generated data
FilePath=fileparts(which(ModEq));
Name=[FilePath,filesep,'data' ModEq(3:end)];
save(Name,'t','datos','proceso')
