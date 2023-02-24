function PlotModelRealization(t,Z,proceso,IND_OBS_X) 
%
% plots the simulated state variable and their observations
%
%         t : time
%         Z : observed data
%   proceso : simulated variables
% IND_OBS_X : Vector of boolean variables defining whether each vector
%             component of the state variable have been observed or not

n_v_e = size(proceso,1);
n_obs = sum(IND_OBS_X);
if n_obs ~= size(Z,1)
    error('Length of vector "IND_OBS_X" wrong')
end

%plot of the state variables
for i=1:n_v_e
   subplot(n_v_e+n_obs,1,i)
   plot(t,proceso(i,:))
   ylabel(['x',num2str(i)])
end

% plot of the observation
for i=1:n_obs
  subplot(n_v_e+n_obs,1,n_v_e+i)
  plot(t,Z(i,:));  
  ylabel(['z',num2str(i)])
end

% figure information
xlabel('t')
subplot(n_v_e+n_obs,1,1); title('generation of data by integrating the SSM')
