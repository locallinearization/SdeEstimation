function [phi, m_new, m, nexpo, err] = phi1_adapt(A,b,h,y,m,rtol,atol,kdmax,kdmin)
% phi1_adapt computes the action of the operator phi_fun(z)=(exp(z)-1)/z over a vector
%
% y + h*phi_fun(h*A)*b = y + phi(:,1)
%
% y + 2*h*phi_fun(2*h*A)*b = y + phi(:,2)
%
% input:
%   A : squared matrix
%   b : vector
%   h : positive number
%   y : reference vector
%   m : initial Krylov dimension
%   rtol: relative tolerance
%   atol: absolute tolerance
%   kdmax: maximun Krylov dimension
%   kdmin: minimun Krylov dimension
%
% output:
%   phi(:,1) = 0.5*h*phi_fun(0.5*h*A)*b
%   phi(:,2) =     h*phi_fun(h*A)*b
%   m_new : Krylov dimension to be used in the next call to k_phi_adapt
%   m : Krylov dimension used to compute phi
%   nexpo : number of computed matrix exponentials
%   err : relative error of the computed phi with respect to the vector y

% Reference
%   "Locally Linearized Runge-Kutta method of Dormand and Prince 
%    for large systems of initial value problems"
%    F.S.Naranjo-Noda, J.C.Jimenez
%    Journal of Computational Physics 426 (2021) 109946 

n=length(b);
m=min(m,n);

% the coeficient is in rfacs( m )
% but is expensive call a function
% and coeficient are hardcode
rfacmin=1;
fac=1/log(2);
gamma = 0.005;

nexpo=0;

tolr= [1.0e-6 1.0e-9 1.0e-12 1.0e-15];
Table=[2 3 4 4;3 4 5 6];

btol = 2*eps;
mb = m;
k1 = 3;
ww = atol+rtol.*abs(y);

V = zeros(n,m+1);
H = zeros(m+3,m+3);

A=h.*A;
rndoff= eps;
normb = norm(b);


% build base of hA and b

% begin Arnoldi
V(:,1)=(1/normb).*b;
% test if norm(b) is 0
if any(isnan(V(:,1)))
    phi = sparse(n,5);
    err = 0;
    warning('llint:phi1LLDP_hJ_f_na',...
       'First vector of krylov subspace seems to be 0.')
    return;
end

for j = 1:m
    p = A*V(:,j);
    s = V(:,1:j);
    H(1:j,j) = s.'*p;
    p = p - s*H(1:j,j);
    s = norm(p);
    if s < btol
        k1 = 0;
        mb = j;
        break
    end
    H(j+1,j) = s;
    V(:,j+1) = (1/s).*p;
end
% end Arnoldi
% using scaling invariance property of Arnoldi
% rescaling H to have H for the subspace of A and b
H=(1/h).*H;
% build \hat{H}
if k1 == 0
    if mb>1
        mb=mb-1;
        warning('llint:phi1LLDP_hJ_f_na',...
       'Breakdown at dimension 1.')
    end
    m=mb;
    hk = H(m+1,m);
    H=[H(1:m,1:m),zeros(m,3);zeros(3,m+3)];
    % V=V(:,1:m+1);
else
    hk = H(m+1,m);
    H(m+1,m) = 0;
end
H(1,m+1) = 1;
H(m+1,m+2) = 1; H(m+2,m+3) = 1;
avm1dot=A*V(:,m+1);

work=1;

while work
    % select p-p of Pade
    nhC = h*norm(H,'inf');
    col = find(tolr>=rtol,1,'last');
    if isempty(col)
      pd = 2;
    else
      fil = (nhC>=1)+1;
      pd = Table(fil,col);
    end
    
    % scaling calculation
    [~,e] = log2(nhC);
    s = max(0,e+1);
    
    % exponential calculation
    M1 = expm64v4(h.*H,pd,s);
    M2 = M1*M1;
    nexpo = nexpo + 1;
    
    % calating \hat{E}
    beta = normb;
    % error relative
    % the divsion by h is because Av_{m+1} is in reality
    % hAv{m+1}(avm1dot) and need rescaling
    err_rel=sqrt((1/n)*sum((((hk*M2(m,m+3)*beta/h).*avm1dot)./ww).^2));
    if err_rel/gamma>=1 && m<kdmax && k1~=0
        % the coeficient is in rfacs( m )
        % but is expensive call a function
        % and coeficient are hardcode
        %  [ rfacmin,rfacmax ] = rfacs( m );
        rfacmax=max(1,m/3);
        knew =  log(err_rel/gamma)*fac;
        knew = ceil(m + min(rfacmax,max(knew,rfacmin)));
        knew = max(kdmin, min(kdmax,knew));

        H = [h.*H(1:m,1:m),zeros(m,knew-m+3);zeros(knew-m+3,knew+3)];
        H(m+1,m) = hk*h;
%         if size(V,2)<knew+1
            V = [V,zeros(n,knew-size(V,2)+1)];
%         end
        mtemp=m+1;
        m=knew;
        mb=m;
        k1 = 3;
        
        j=mtemp;
        s = V(:,1:j);
        H(1:j,j) = s.'*avm1dot;
        avm1dot = avm1dot - s*H(1:j,j);
        s = norm(avm1dot);
        if s < btol
            k1 = 0;
            mb = j;
            mtemp=m;
        end
        H(j+1,j) = s;
        V(:,j+1) = (1/s)*avm1dot;

        for j = mtemp+1:m
            p = A*V(:,j);
            s = V(:,1:j);
            H(1:j,j) = s.'*p;
            p = p - s*H(1:j,j);
            s = norm(p);
            if s < btol
                k1 = 0;
                mb = j;
                break;
            end
            H(j+1,j) = s;
            V(:,j+1) = (1/s)*p;
        end
        H=(1/h).*H;
        % build \hat{H}
        if k1 == 0
            if mb>1
                mb=mb-1;
            end
            m=mb;
            hk = H(m+1,m);
            H=[H(1:m,1:m),zeros(m,3);zeros(3,m+3)];
        else
            hk = H(m+1,m);
            H(m+1,m) = 0;
        end
        H(1,m+1) = 1;
        H(m+1,m+2) = 1; H(m+2,m+3) = 1;
        avm1dot=A*V(:,m+1);
    else
        rfacmax = - m/4;
        rfacmin = m/3;
        m_new =  log(err_rel/gamma)*fac;
        m_new = floor(m + max(rfacmax,min(m_new,rfacmin)));
        m_new = max(kdmin, min(kdmax,m_new));
        break
    end
end

mx = m + 1;

M1(mx,mx) = hk*M1(m,m+2);
M1(m+2,mx) = hk*M1(m,m+3);
M2(mx,mx) = hk*M2(m,m+2);
M2(m+2,mx) = hk*M2(m,m+3);

avm1dot = V(:,1:mx);
% Matix projection
phi = avm1dot*(beta*[M1(1:mx,mx) M2(1:mx,mx)]);

err = max(err_rel,rndoff);
end

