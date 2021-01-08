function [x,u,err,Psi,Phi]=scalar_hyp(xa,xb,t0,T,beta,uex,u0,ul,nx,deltat,param)
% SCALAR_HYP Numerical solution of scalar linear hyperbolic equations,
%     formulas (3.7.1a), pag. 145, CHQZ2
%
%      du/dt +beta du/dx=0           in (xa,xb) x (0,T). 
%      initial cond.: u(x,t)=u0(x),  
%      inflow cond.: u(xa,t)=ul(t)
%
%  Space discretization: 1. Galerkin-Numerical Integration with LGL formulas,
%                        2. Legendre collocation (strong form),
%                        3. penalty,
%
%  Time discretization:  RK4
%
%  [x,u,err,Psi,Phi]=scalar_hyp(xa,xb,t0,T,beta,uex,u0,ul,nx,deltat,param);
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
%        param(1)= space disc. scheme: 1=strong (3.7.2), pag. 146 CHQZ2, 
%                                     2=GNI+ weak b.c. (3.7.5), pag. 147
%                                     3=penalty (3.7.2a)+(3.7.9), pag. 148
%        param(2)=tau : penalty parameter (used only if param(1)==3)
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
% param(1)=2;param(2)=0.75
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
jac=(xb-xa)*5.d-1;
x=x*jac+(xa+xb)*5.d-1;
alpha=25/12; alpha0=4; alpha1=-3; alpha2=4/3; alpha3=-0.25;


tau=param(2)/(jac*wx(1));
t=t0;
un=zeros(npdx,1); un=u0(x);
intu0=sum(un.*wx)*jac;
UL=[UL;un(npdx)-ul(t)];
INTU=[INTU;0];
intu1ul=0;
psi=0;
Psi=[Psi;psi];




% Construction of matrix A

if param(1)==1
% strong collocation (3.7.2)
A=-beta*dx/jac;

elseif param(1)==2
% GNI+weak b.c. (3.7.5)
A=beta*dx'*diag(wx);
A(npdx,npdx)=A(npdx,npdx)-beta;

elseif param(1)==3
% penalty (3.7.2a) +(3.7.9)
A=-beta*diag(wx)*dx;
A(1,1)=A(1,1)-tau*beta*wx(1)*jac;

end
deltat2=deltat*0.5;




%close all
% fig=figure(...,
%     'Name','Hyperbolic 1d problem',...
%     'Visible','on')
% axis([-1,1,-1,1]); 
% time loop explicit 4th order Runge Kutta 
for it=1:nt-1
t=tt(it);
if(param(1)==1)
w1=A*un;
w2=A*(un+deltat2*w1);
w3=A*(un+deltat2*w2);
w4=A*(un+deltat*w3);
u=un+deltat/6*(w1+2*w2+2*w3+w4); u(1)=ul(tt(it+1));
elseif(param(1)==2)
w1=A*un; w1(1)=w1(1)+beta*ul(t); w1=(w1./wx)/jac;
w2=A*(un+deltat2*w1); t=t+deltat2; w2(1)=w2(1)+beta*ul(t);w2=(w2./wx)/jac;
w3=A*(un+deltat2*w2);w3(1)=w3(1)+beta*ul(t);w3=(w3./wx)/jac;
w4=A*(un+deltat*w3); t=t+deltat2; w4(1)=w4(1)+beta*ul(t);w4=(w4./wx)/jac;
u=un+deltat/6*(w1+2*w2+2*w3+w4);
elseif(param(1)==3)
w1=A*un; w1(1)=w1(1)+tau*beta*ul(t)*wx(1)*jac;w1=(w1./wx)/jac;
w2=A*(un+deltat2*w1); t=t+deltat2; w2(1)=w2(1)+tau*beta*ul(t)*wx(1)*jac;
w2=(w2./wx)/jac;
w3=A*(un+deltat2*w2);w3(1)=w3(1)+tau*beta*ul(t)*wx(1)*jac;w3=(w3./wx)/jac;
w4=A*(un+deltat*w3); t=t+deltat2; w4(1)=w4(1)+tau*beta*ul(t)*wx(1)*jac;
w4=(w4./wx)/jac;
u=un+deltat/6*(w1+2*w2+2*w3+w4);
end

intu=sum(u.*wx)*jac;
t=tt(it+1); UL=[UL;u(npdx)-ul(t)];
INTU=[INTU;intu-intu0];

if(fix((it+1)/2)*2~=it+1)
intu1ul=intu1ul+simpcz(tt(it-1:it+1),UL(it-1:it+1));
psi=intu-intu0+beta*intu1ul;
Psi=[Psi;psi];
end
nor=sum(u.*u.*wx)*jac;

% plot(x,u); 
% title(['t=',num2str(t)]);
% pause(0.01)
% 

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
u1=u(1);u2=u(npdx);
phi1=5.d-1*beta*(-u1^2+u2^2)+...
5.d-1*delt1*(alpha*nor-(alpha0*nor0+alpha1*nor1+alpha2*nor2+alpha3*nor3));
Phi=[Phi;phi1];
un3=un2; un2=un1; un1=un; un=u;
nor3=nor2;nor2=nor1;nor1=nor0; nor0=nor;
end
end


% error  at  time T

ue=uex(x,tt(nt));
err=norm(u-ue);

% fig=figure(...,
%     'Name','Hyperbolic 1d problem Phi(t)',...
%     'Visible','on')
% plot(tt(5:nt-1),Phi);
% ylabel('Phi(t)'); xlabel('t')
