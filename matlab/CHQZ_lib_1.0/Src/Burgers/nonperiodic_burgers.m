function [u,err]=nonperiodic_burgers(xa,xb,t0,T,visc,nx,deltat,tscheme,bc)
% NONPERIODIC_BURGERS Numerical solution of non periodic Burgers equation
%      
%  \partial u/\partial t +u \partial u/\partial x -\nu \partial^2 u/\partial x^2
%                 =0     in Omega, \forall t>t_0
%  u(x,t_0)=u_0(x)    in Omega
%  + b.c.
%
%     Exact solution is given by formula (3.1.8), pag. 119, CHQZ2
%     Space discretization: Galerkin-Numerical Integration with LGL
%     formulas.
%     Time discretization:  Explicit Euler, RK2, RK4
%
%  [u,err]=nonperiodic_burgers(xa,xb,t0,T,visc,nx,deltat,tscheme,bc)
%
% Input: xa,xb  = extrema of space domain Omega=(xa,xb)
%        t0,T   = extrema of time domain [t0,T]
%        visc   = viscosity
%        nx     = spectral polynomial degree, the number of LGL points is
%                 nx+1
%        deltat = time step
%        tscheme= time discretization scheme: 1=EE, 2=RK2, 3=RK4
%        bc     = choice of boundary conditions: 1 == Dirichlet 
%                                                2 == Neumann 
%
% Output: u = solution at final time T
%         err= L_\infty norm error between numerical and exact solution at
%         final time T
%
%
% Possible input data:
% xa=-1;xb=1;
% visc=0.01;
% tscheme=4;bc=1;
% nx=16;
% t0=0;T=1;
% deltat=1.e-4;
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$


%  Definition of exact solution and first derivative to set Neumann b.c
 
Uex=['1+(x-t)/(t+1).*(4/sqrt(t+1)*exp(-(x-t).^2/(0.04*(t+1))))',...
    './(1+4/sqrt(t+1)*exp(-(x-t).^2/(0.04*(t+1))))'];
Uex1=['-4*(50*(t+1)^(1/2)*t^2-100*x*t*(t+1)^(1/2)-4*exp(-25*(x-t)^2/(t+1))*t',...
'+50*x^2*(t+1)^(1/2)-(t+1)^(1/2)*t-(t+1)^(1/2)-4*exp(-25*(x-t)^2/(t+1)))',...
    '*exp(-25*(x-t)^2/(t+1))/((t+1)^(1/2)+4*exp(-25*(x-t)^2/(t+1)))^2/(t+1)^2'];
uex=vectorize(Uex);
uex1=vectorize(Uex1);
uex=@(x,t)[eval(uex)];
uex1=@(x,t)[eval(uex1)];



% space discretization

npdx=nx+1; 
[x,w]=xwlgl(npdx); 
dx=derlgl(x,npdx);
jac=(xb-xa)*5.d-1;
x=x*jac+(xb+xa)*5.d-1;
A=stiff_1d_sp(w,dx,jac);
A=visc*A;

t=t0;
u0=uex(x,t);

fig=figure(...
    'Name','Non periodic Burgers solution',...
    'Visible','on');
title(['t=',num2str(t)]);
%U1=uex(x);
%plot(x,u0,x,U1,'r',x,zeros(npdx,1),'+');
%pause(0.001)
k=0;

% time loop

while t< T
if tscheme==1
    % Explicit Euler in time
u=ee(t,deltat,@fburgers,u0,dx,A,w,visc,uex,uex1,bc);
elseif tscheme==2
    % RK2 in time
u=rk2(t,deltat,@fburgers,u0,dx,A,w,visc,uex,uex1,bc);
elseif tscheme==4
    % RK4 in time
u=rk4(t,deltat,@fburgers,u0,dx,A,w,visc,uex,uex1,bc);
end
t=t+deltat;
k=k+1;
U1=uex(x,t);
% if rem(k,100)==0
% err=norm(u(1:nx)-U1(1:nx),inf);
% fprintf('t=%13.6e, nx=%d, err_inf=%13.6e, tscheme=%d \n ',t,nx, err,tscheme)
% end
plot(x,u,x,U1,'r',x,ones(npdx,1),'+');
axis([-1,1,0.8,1.2])
title(['t=',num2str(t)]);
pause(0.001)
u0=u; 
if (sum(isnan(u))>1)
    fprintf('Nan solution: instability due to too large time-step')
    err=100;
    return
end
% end time loop
end
U1=uex(x,t);
err=norm(u(1:nx)-U1(1:nx),inf);
fprintf('t=%13.6e, nx=%d, err_inf=%13.6e, tscheme=%d \n ',t,nx, err,tscheme)


%U1=uex(x,t);
%err=norm(u(1:nx)-U1(1:nx),inf);
%disp('N, error')
%[nx, err]
plot(x,u,x,U1,'r',x,ones(npdx,1),'+');
axis([-1,1,0.8,1.2])
title(['t=',num2str(t)]);
%pause(0.001)

return
