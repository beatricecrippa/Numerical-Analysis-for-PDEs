function [f]=fburgers(tt,deltat,u0,A,dx,w,visc,uex,uex1,bc);
% FBURGERS Defines the non-per Burgers function for nonperiodic_burgers.m
%
% [f]=fburgers(tt,deltat,u0,A,dx,w,visc,uex,uex1,bc);
%
% Input : tt =time
%         deltat = time step
%         u0 = initial data
%         A = matrix related to differential operator
%         dx = LGL first derivative matrix
%         w = LGL weights array
%         visc = viscosity coefficient (constant > 0)
%         uex  = exact solution (uex=@(x)[uex(x)], with .*, .^, ./)
%         uex1 = exact solution (uex=@(x)[uex(x)], with .*, .^, ./)
%         bc     = choice of boundary conditions: 1 == Dirichlet 
%                                                 2 == Neumann 
%
% Output: f = array with evaluation of fburgers
%   Written by Paola Gervasio
%   $Date: 2007/04/01$
%

np=length(u0);
b=zeros(np,1);
if bc ==2
t=tt;
uu1=uex(-1); ux1=uex1(-1);
b(1)=0.5*uu1^2-visc*ux1;
uu1=uex(1); ux1=uex1(1);
b(np)=-0.5*uu1^2+visc*ux1;
end
f=(b-A*u0+0.5*dx'*diag(w)*(u0.*u0))./w;

return

