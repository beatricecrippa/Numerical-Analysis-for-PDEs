function [u1]=rk2(t,deltat,f,u0,dx,A,w,visc,uex,uex1,bc);
% RK2 One step of  Explicit Runge-Kutta second order scheme.
% 
% [u1]=rk2(tt,deltat,f,u0,dx,A,w,visc,uex,uex1,bc);
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
% Output: u1 = array solution, computed with one step of RK2
%

%   Written by Paola Gervasio
%   $Date: 2007/04/01$
%

K1=feval(f,t,deltat,u0,A,dx,w,visc,uex,uex1);
t1=t+deltat; u1=u0+deltat*K1;
if bc ==1
u1(1)=uex(-1,t1); u1(np)=uex(1,t1);
end

K2=feval(f,t1,deltat,u0,A,dx,w,visc,uex,uex1);
u1=u0+deltat/2*(K1+K2);
if bc ==1
u1(1)=uex(-1,t1); u1(np)=uex(1,t1);
end

return

