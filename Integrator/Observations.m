function z = Observations(t,x,theta,FileNames,ObsX,seed)
%
% Generates a time series z of noisy observations of the state varible x of the model
%     dx = f(t,x,theta) dt + g(t,x,theta) dw, 
% with the discrete observation equation 
%      z = h(t,x,theta) + p(t,x,theta)n + e.
%
% Notation
%   n_m: number of time-instants where "x" is observed
%   n_obs: number of observation equations
%   n_v_e: number of state variables "x" 
%   n_theta: number of parameters theta 
%
%  Imput
%             t: integration times                                                 (1 x n_m)
%             x: numerical solution of the state variable                          (n_v_e x n_m)
%         theta: model parameters                                                  (1 x n_theta)
%     FileNames: cell variable with the name of the functions defining the model   (1 x 3)
%          ObsX: index of the observed state variables                             (n_v_e x 1)
%          seed: seed for initializing the random number generator                 (1 x 1)
%
%  Output
%             z: noisy observations                                                (n_obs x n_m)

ObsEq = FileNames{2};

if (nargin==6) & (~isempty(seed))
   randn('seed',seed);
end

n_obs = sum(ObsX);
[~,n_m] = size(x);

if nargout(ObsEq)==3 
   z=zeros(n_obs,n_m);  
   for k=1:n_m
     [he,~,See] = ObsEq(t(k),x(:,k),theta);
     if norm(See)<eps
       SeeRc=See;
     else   
       SeeRc = sqrtm(See);
     end   
     z(:,k) = he + SeeRc * randn(n_obs,1);
   end
end

if nargout(ObsEq)==6
   z=zeros(n_obs,n_m); 
   for k=1:n_m
     [he,~,See,pe,~,Lee] = ObsEq(t(k),x(:,k),theta);
     if norm(See)<eps
       SeeRc=See;
     else   
       SeeRc = sqrtm(See);
     end   
     n=size(pe);
     for i=1:n
       if norm(Lee(i))<eps
         LRc(i)=Lee(i);
       else   
         LRc(i) = sqrt(Lee(i));
       end   
     end
     z(:,k) = he + SeeRc * randn(n_obs,1);
     for i=1:n
        z(:,k)= z(:,k) + pe{i} * LRc(i) * randn;
     end  
   end  
end

