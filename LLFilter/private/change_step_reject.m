function hnew=change_step_reject(h,error,pp)
facmax=1;
facmin=0.1;
%fac=0.25;
fac=0.2;
hnew=h.*min(facmax,max(facmin,fac.*(1./error).^pp));
