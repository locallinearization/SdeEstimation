function [Xf,Pf,Xp,Pp,Inn,PInn,pts,LogLK]=LLfilter_Det(t,Z,Xf0,Pf0,theta,FileNames,LKWindow,max_pts,AbsTol,RelTol)
%    Adaptive Deterministic Local Linearization Filter for the estimation of the 
% unobserved componet x of the continuous model
%        dx = f(t,x,theta) dt + g(t,x,theta) dw, 
% with dicrete observation equation 
%         Z = h(t,x,theta) + p(t,x,theta)n + e, 
% given a set of discrete observation of Z on t. 
% Xf0, Pf0 and theta are also known and given in advance.
%
% Notation
%   n_m: number of time-instants where Z is observed
%   n_obs: number of observation equations
%   n_v_e: number of state variables x 
%   n_theta: number of parameters theta 
%
% Input:
%   t: instants of times where at least one component of Z is observed  (1 x n_m)          
%   Z: observed sample                        (n_obs x n_m)   Set Z(k,t)=inf when data Z(k,t) is missing
%   Xf0: initial filter value                 (n_v_e x 1)
%   Pf0: initial filter variance value        (n_v_e x n_v_e)
%   theta: State Space Model parameters       (1 x n_theta)
%   FileNames={ModEq,ObsEq,ExactSol}   
%      ModEq string with the .m file name of the continuous model definition
%      ObsEq string with the .m file name of the observation equation definition
%      ExactSol string with the .m file name of the exact E(x/Zn) and E(xx'/Zn) 
%   LKWindow: interval of time where the likelihood function will be computed  (1 x 2)
%             set LKWindow=[] to get LKWindow = [t(1) t(end)]
%   max_pts: maximum number of allowed misssing points to compute the
%            prediction between two observations   (1 x n_m)
%            max_pts=[] implies not restriction
%   AbsTol: Absolute Tolerance (1 x 3), 1: for the first moment, 
%                                       2: for the second moment, 
%                                       3: for the number of simulations
%   RelTol: Relative Tolerance (1 x 3), 1,2,3 same than AbsTol 
%
% OutPuts:
%   Xf : LL filter  
%   Pf : LL filter variance 
%   Xp : LL prediction 
%   Pp : LL prediction variance  
%   Inn: Innovation
%   PInn: Innovation variance
%   pts: number of misssing points used to compute the prediction between two observations (2 x n_m)
%         row 1: accepted, row 2: failed
%   LogLK:  -2log(Likelihood)
%

%  References
%   - "Approximate linear minimum variance filters for continuous-discrete 
%      state space models: convergence and practical adaptive algorithms", 
%      Jimenez J.C.
%      IMA Journal of Mathematical Control and Information, 
%      Volume 36, Issue 2, June 2019, Pages 341–378. 
%   - "Simplified formulas for the mean and variance of linear 
%      stochastic differential equations",  
%      Jimenez J.C.
%      Applied Mathematics Letters
%      Volume 49, November 2015, Pages 12-19. 

if ~isempty(max_pts) && (length(max_pts)~=length(t))
    error('max_pts wrong. Check the lenght of max_pts')
end

n_m=size(t,2);        % number of instants of times
n_v_e=size(Xf0,1);    % number of state variables
nobs=size(Z,1);       % number of observed variables
if nargout(FileNames{1})==2  % pp=1/(1+order)   with order= 1 or 2 
    pp=1/3; 
else
    pp=0.5;
end
rtol=min(RelTol(1:2));
atol=min(AbsTol(1:2));

Xf=zeros(n_v_e,n_m);
Xp=zeros(n_v_e,n_m);
Pf=zeros(n_v_e,n_v_e,n_m);
Pp=zeros(n_v_e,n_v_e,n_m);
Inn=[zeros(nobs,1) inf*ones(nobs,n_m-1)];
PInn=zeros(nobs,nobs,n_m);
pts=zeros(2,n_m);

Zt=Xf0;
Vt=Pf0;
Xf(:,1)=Zt;
Xp(:,1)=Zt;
Pf(:,:,1) = Pf0;
Pp(:,:,1) = Pf0;

LogLK = 0.0;
[Prediction, Filter, ModEq, ObsEq] = ModObsEqType(FileNames,n_v_e);
hh=initial_step(t(1),Zt,Vt + Zt*Zt',theta,ModEq,AbsTol,RelTol,pp);
max_pts = max_pts - 1;
for k= 2:n_m
   
  Ztdt=Zt;
  Pt  =Vt + Zt*Zt';
  
  ti=t(k-1);
%  hh=initial_step(ti,Ztdt,Pt,theta,FileNames,AbsTol,RelTol,pp);
  
  follow = 1;
  while (t(k)-ti>eps) && follow
     TT=0.5.*(t(k)-ti);
     if ~isempty(max_pts) && ((pts(1,k)+pts(2,k)) == max_pts(k)), hh=TT; follow=0; end  % to scape of nonconvergen realization
     if hh>TT, hh=TT; end
     hh2=2.*hh;
     ti_hh =ti+hh;
     ti_hh2=ti+hh2;
     [Ztdt1, Pt1, Ztdt2, Pt2] = Prediction(ti_hh,ti,Ztdt,Pt,theta,ModEq,rtol,atol);  
     [Ztdt3, Pt3] = Prediction(ti_hh2,ti_hh,Ztdt1,Pt1,theta,ModEq,rtol,atol);  
     errorZ=error_step(Ztdt,Ztdt3,Ztdt2,n_v_e,AbsTol(1),RelTol(1));
     errorP=error_step(vec(Pt),vec(Pt3),vec(Pt2),sqr(n_v_e),sqr(AbsTol(2)),RelTol(2));
     if (errorZ<=1) && (errorP<=1)
        ti=ti_hh2;
        Ztdt=Ztdt3;
        Pt=Pt3;
        hhZ=change_step_no_reject(hh,errorZ,pp);
        hhP=change_step_no_reject(hh,errorP,pp);
        pts(1,k)=pts(1,k)+1;
     else
        if errorZ>1
           hhZ=change_step_reject(hh,errorZ,pp);
        else
           hhZ=change_step_no_reject(hh,errorZ,pp);         
        end
        if errorP>1
           hhP=change_step_reject(hh,errorP,pp);
        else   
           hhP=change_step_no_reject(hh,errorP,pp);
        end
        pts(2,k)=pts(2,k)+1; 
     end   
     hh=max(eps,min([hhZ hhP]));
  end      
  Ut = Pt - Ztdt*Ztdt';
  
  [Zt, Vt, nu, Sv] = Filter(t(k),Ztdt,Ut,theta,Z(:,k),ObsEq);    % LL filter & innovation 
  
  if (LKWindow(1)<=t(k)) && (t(k)<=LKWindow(2))
     LogLK = LogLK + log(det(Sv)) + nu' /Sv * nu;                  
  end

  Xp(:,k)   = Ztdt;    % LL prediction
  Pp(:,:,k) = Ut;      % LL prediction variance
  Xf(:,k)   = Zt;      % LL filter 
  Pf(:,:,k) = Vt;      % LL filter variance
  
  obs= find(Z(:,k)~=inf);
  Inn(obs,k)  = nu;      % innovation
  PInn(obs,obs,k) = Sv;  % innovation variance
   
end

LogLK = LogLK + log(2.*pi).*n_m;


%----------------------------

function v = vec(x)

[n, m]=size(x);
v=reshape(x,n*m,1);


%----------------------------
% LL filter initial step size 

function  hh=initial_step(t0,m0,P0,theta,FileName,AbsTol,RelTol,pp)

d=size(m0,1);
switch nargout(FileName)
    case 2
        [f0,A] = FileName(t0,m0,theta);
        a1=zeros(d,1);
        g0{1}=a1;
    case 3
        [f0,g0,A] = FileName(t0,m0,theta);
        a1=zeros(d,1);
    case 4
        [f0,g0,A,B] = FileName(t0,m0,theta);
        a1=zeros(d,1);
    case 5
        [f0,g0,A,a1,b1] = FileName(t0,m0,theta);
    case 6
        [f0,g0,A,B,a1,b1] = FileName(t0,m0,theta);
    otherwise
        disp('Error in Model definition')
end

hh_m=starting_step(f0,A,a1,m0,d,AbsTol(1),RelTol(1),pp);

d2=d*d;
I_d=eye(d);
bbT=zeros(d,d);
b1bT=zeros(d,d);
bb1T=zeros(d,d);
BK{1}=zeros(d2,d2);
BK{2}=zeros(d2,d);
BK{3}=zeros(d2,d);
BK{4}=zeros(d2,d);
BK{5}=zeros(d2,d);
m=size(g0,2);
if nargout(FileName)==3 || nargout(FileName)==5 || nargout(FileName)==2    %B=0
   if nargout(FileName)==3 || nargout(FileName)==2 %b1=0;
      for i=1:m 
        bbT = bbT + g0{i}*g0{i}';
      end
   else
      for i=1:m 
        bbT = bbT + g0{i}*g0{i}';
        b1bT = b1bT + b1{i}*g0{i}';
        bb1T = bb1T + g0{i}*b1{i}';
      end
   end
else
   b  = bxt(m0,g0,B);
   if nargout(FileName)==4 %b1=0;
     for i=1:m 
       bbT = bbT + b{i}*b{i}';
       BK{1}= BK{1} + kron(B{i},B{i});
       BK{2}= BK{2} + kron(b{i},B{i});
       BK{3}= BK{3} + kron(B{i},b{i});
     end
   else
     for i=1:m 
       bbT = bbT + b{i}*b{i}';
       BK{1}= BK{1} + kron(B{i},B{i});
       BK{2}= BK{2} + kron(b{i},B{i});
       BK{3}= BK{3} + kron(B{i},b{i});
   
       BK{4}= BK{4} + kron(b1{i},B{i});
       BK{5}= BK{5} + kron(B{i},b1{i});
       b1bT = b1bT + b1{i}*b{i}';
       bb1T = bb1T + b{i}*b1{i}';
     end
   end
end
a0=f0-A*m0;
AE = kron(I_d,A) + kron(A,I_d) + BK{1};
beta1 = vec(bbT);
beta2 = vec(b1bT + bb1T);
beta4 = kron(I_d,a0) + kron(a0,I_d) + BK{2} + BK{3};
beta5 = kron(I_d,a1) + kron(a1,I_d) + BK{4} + BK{5};
F0=AE*vec(P0)+beta1+beta4*m0;
Ft0=beta2+2*beta5*m0+beta4*f0;
hh_P=starting_step(F0,AE,Ft0,vec(P0),d2,AbsTol(2),RelTol(2),pp);

hh = max(eps,min([hh_m; hh_P]));


