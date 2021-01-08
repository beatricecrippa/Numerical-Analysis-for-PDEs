function [x,u,err,Psi,Phi]=stag_scalar_hyp(xa,xb,t0,T,beta,uex,u0,ul,nx,deltat)
% STAG_SCALAR_HYP  Numerical solution of scalar linear hyperbolic equations,
%     formulas (3.7.1a), pag. 145, CHQZ2 by staggered grids method (pag. 149, 
%     CHQZ2)
%
%      du/dt +beta du/dx=0           in (xa,xb) x (0,T). 
%      initial cond.: u(x,t)=u0(x),  
%      inflow cond.: u(xa,t)=ul(t)
%
%  Space discretization: Staggered grids method (pag. 149, CHQZ2)
%
%  Time discretization:  RK4
%
%  [x,u,err,Psi,Phi]=stag_scalar_hyp(xa,xb,t0,T,beta,uex,u0,ul,nx,deltat,param);
%
% Input: xa,xb  = extrema of space domain Omega=(xa,xb)
%        t0,T   = extrema of time domain [t0,T]
%        beta   = transport coefficient
%        uex    = exact solution (uex=@(x)[uex(x)], with .*, .^, ./)
%        u0    = initial condition (u0=@(t)[u0(t)], with .*, .^, ./)
%        ul    = inflow condition (ul=@(x)[ul(x)], with .*, .^, ./)
%        nx     = spectral polynomial degree, the number of LGL points is
%                 nx+1
%        deltat = time step
%
% Output: u = solution at final time T
%         err= L_\infty norm error between numerical and exact solution at
%                    final time T
%         Psi = column array of evaluations Psi(t_n) (see def. 3.7.14),
%               pag. 151, CHQZ2
%         Phi = columen array of evaluations Phi(t_n) where
%               Phi(t)=1/2 d/dt ||u||^2+1/2beta*(u^2(1)-u^2(-1))
%
% Possible input data:
% xa=-1;xb=1;
% beta=1.5;
% nx=8;
% t0=0;T=4;
% deltat=1.e-4;
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$


delt1=1/deltat;
tt=(t0:deltat:T)';nt=length(tt);
Phi=[]; UN1=[]; UL=[];Psi=[]; INTU=[];
npdx=nx+1;
[x,wx]=xwlgl(npdx); 
[dx]=derlgl(x,npdx);
[xg,wxg]=xwlg(nx);
jac=(xb-xa)*5.d-1;
alpha=25/12; alpha0=4; alpha1=-3; alpha2=4/3; alpha3=-0.25;

% eta = interpolation matrix from LG grid (xg) to LGL grid (x)

[eta]=intlag_lg(xg, wxg, x);

% ltrm = interpolation matrix from LGL grid to LG grid

[ltrm]=intlag_lgl(x, xg);

% dg Legendre Gauss derivative matrix 

[dg] = derlg(xg,nx); dg=dg/jac;


x=x*jac+(xa+xb)*5.d-1;
xg=xg*jac+(xa+xb)*5.d-1;


t=t0;
un=zeros(nx,1); un=u0(xg);
intu0=sum(un.*wxg)*jac;u2=eta(npdx,:)*un;
UL=[UL;u2-ul(t)];INTU=[INTU;0];
intu1ul=0;
psi=0;
Psi=[Psi;psi];


utilde=zeros(npdx,1);


deltat2=deltat*0.5;

for it=1:nt-1
t=tt(it);
% 
utilde=eta*un;
 utilde(1)=ul(t);
%
w1=-beta*dg*(ltrm*utilde);

% 
utemp=un+deltat2*w1;
utilde=eta*utemp;
t=t+deltat2; utilde(1)=ul(t);
w2=-beta*dg*(ltrm*utilde);

% 
utemp=un+deltat2*w2;
utilde=eta*utemp;
utilde(1)=ul(t);
w3=-beta*dg*(ltrm*utilde);

% 
utemp=un+deltat*w3;
utilde=eta*utemp;
t=t+deltat2; utilde(1)=ul(t);
w4=-beta*dg*(ltrm*utilde);

u=un+deltat/6*(w1+2*w2+2*w3+w4);
nor=sum(u.*u.*wxg)*jac;
t=tt(it+1); u1=ul(t);u2=eta(npdx,:)*u;
intu=sum(u.*wxg)*jac;
UL=[UL;u2-u1];INTU=[INTU;intu-intu0];
if(fix((it+1)/2)*2~=it+1)
intu1ul=intu1ul+simpcz(tt(it-1:it+1),UL(it-1:it+1));
psi=intu-intu0+beta*intu1ul;
Psi=[Psi;psi];
end


% evaluate Phi(t)=1/2 d/dt ||u||^2+1/2beta*(u^2(1)-u^2(-1))
if(it==1)
un=u; nor0=nor;
elseif (it==2)
un1=un; un=u;nor1=nor0; nor0=nor;
elseif(it==3)
un2=un1; un1=un; un=u;nor2=nor1;nor1=nor0; nor0=nor;
elseif(it==4)
un3=un2; un2=un1; un1=un; un=u;
nor3=nor2;nor2=nor1;nor1=nor0; nor0=nor;
else
phi1=5.d-1*beta*(-u1^2+u2^2)+...
5.d-1*delt1*(alpha*nor-(alpha0*nor0+alpha1*nor1+alpha2*nor2+alpha3*nor3));
Phi=[Phi;phi1];
un3=un2; un2=un1; un1=un; un=u;
nor3=nor2;nor2=nor1;nor1=nor0; nor0=nor;
end
end

% compute error at T
ue=uex(xg,tt(nt));
err=norm(u-ue);

return
