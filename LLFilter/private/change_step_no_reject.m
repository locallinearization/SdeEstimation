function hnew=change_step_no_reject(h,error,pp)
facmax=5;
facmin=0.25;
fac=0.8;
% para ecuaciones lineales es mejor
% facmax=8;
% facmin=0.2;
% fac=0.9;
hnew=h.*min(facmax,max(facmin,fac.*(1./error).^pp));
