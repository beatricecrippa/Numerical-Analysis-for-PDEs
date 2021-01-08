function [u1]=ee(tt,deltat,f,u0,dx,A,w,visc,uex,uex1,bc);
% EE One step of  Explicit Euler scheme for 1D parabolic problems. 
% 
% [u1]=ee(tt,deltat,f,u0,dx,A,w,visc,uex,uex1,bc);
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
% Output: u1 = array solution, computed with one step of EE
%
%   Written by Paola Gervasio
%   $Date: 2007/04/01$
%

K1=feval(f,tt,deltat,u0,A,dx,w,visc,uex,uex1,bc);
u1=u0+deltat*K1;
if bc==1
np=length(u1);
u1(1)=uex(-1,tt);
u1(np)=uex(1,tt);
end

return

