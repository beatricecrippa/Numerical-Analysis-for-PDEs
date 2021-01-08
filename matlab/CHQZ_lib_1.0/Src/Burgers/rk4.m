function [u1]=rk4(tt,deltat,f,u0,dx,A,w,visc,uex,uex1,cb);
% RK4 One step of  Explicit Runge-Kutta fourth order scheme.
% 
% [u1]=rk4(tt,deltat,f,u0,dx,A,w,visc,uex,uex1,bc);
%
% Input : tt =time
%         deltat = time step
%         f =  function for the r.h.s of the parabolic equation
%         u0 = initial data
%         dx = LGL first derivative matrix
%         A = matrix related to differential operator
%         w = LGL weights array
%         visc = viscosity coefficient (constant > 0)
%         uex  = exact solution (uex=@(x)[uex(x)], with .*, .^, ./)
%         uex1 = exact solution (uex=@(x)[uex(x)], with .*, .^, ./)
%         bc     = choice of boundary conditions: 1 == Dirichlet 
%                                                 2 == Neumann 
%
% Output: u1 = array solution, computed with one step of RK4
%

%   Written by Paola Gervasio
%   $Date: 2007/04/01$
%

np=length(u0);
del2=deltat*0.5;
K1=feval(f,tt,deltat,u0,A,dx,w,visc,uex,uex1,cb);
t1=tt+del2; u1=u0+del2*K1;
if cb ==1
u1(1)=uex(-1,t1); u1(np)=uex(1,t1);
end

K2=feval(f,t1,deltat,u1,A,dx,w,visc,uex,uex1,cb);
u1=u0+del2*K2;
if cb ==1
u1(1)=uex(-1,t1); u1(np)=uex(1,t1);
end
K3=feval(f,t1,deltat,u1,A,dx,w,visc,uex,uex1,cb);
t1=tt+deltat; u1=u0+deltat*K3;
if cb ==1
u1(1)=uex(-1,t1); u1(np)=uex(1,t1);
end
K4=feval(f,t1,deltat,u1,A,dx,w,visc,uex,uex1,cb);
u1=u0+deltat/6*(K1+2*K2+2*K3+K4);
if cb ==1
u1(1)=uex(-1,t1); u1(np)=uex(1,t1);
end

return

