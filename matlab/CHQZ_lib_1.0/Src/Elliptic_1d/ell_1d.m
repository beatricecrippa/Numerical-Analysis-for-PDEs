function [xy,un]=ell_1d(xa,xb,nu,beta,gam,ff,cb,ub,ne,nx);
% ELL_1D   Numerical solution of the 1D elliptic boundary value problem (4.3.17) pag. 208, CHQZ2
%
%  -nu d^2 u/dx^2 +beta du/dx+ gam u=1  in (-1,1)
%
%   + Dirichlet or Neumann boundary conditions
%
% by Galerkin Numerical Integration with LGL quadrature formulas.
%
%  [xy,un]=ell_1d(xa,xb,nu,beta,gam,ff,ub,ne,nx);
%
% Input: xa, xb = extrema of computational domain Omega=(xa,xb)
%      nu   = viscosity (constant>0)
%      beta = coefficient of first order term (constant)
%      gam = coefficient of zeroth order term (constant>0)
%      cb = string containing boundary conditions, for example
%           cb='dn' imposes Dirichlet in xa and Neumann in xb,
%           possible b.c.: Dirichlet, Neumann
%      ub = 2 components array with Dirichlet, Neumann (pure derivative)
%           data in xa and in xb
%      ne = number of elements (equally spaced)
%      nx = polynomial degree in each element (the same in each element)
%
% Output: xy = 1D mesh
%         un = numerical solution 
%
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$


% 
npdx=nx+1; 

[x,wx]=xwlgl(npdx); dx=derlgl(x,npdx);

% nov
ldnov=npdx; nov=zeros(ldnov,ne);
[nov]=cosnov_1d(npdx,ne,nov);
noe=nov(npdx,ne);

% Uniform decomposition of  Omega in ne elements

[xx,jacx,xy,ww]=mesh_1d(xa,xb,ne,npdx,nov,x,wx);

% Variable coefficients evaluation


% Matrix assembling
A=ell_1d_se(npdx,ne,nov,nu,beta,gam,wx,dx,jacx);

% Right hand side
f=ff(xy).*ww;


% Setting Dirichlet boundary conditions on the matrix

if cb(1)=='d'
f(1)=ub(1); i1=2;
elseif cb(1)=='n'
f(1)=f(1)-nu*ub(1);i1=1;
end
if cb(2)=='d';
f(noe)=ub(2);i2=noe-1;
elseif cb(1)=='n'
f(noe)=f(noe)+nu*ub(2);i2=noe;
end

% Reduce the system to non-Dirichlet unknowns
lint=(i1:i2);
if cb=='dd'
f(lint)=f(lint)-A(lint,1)*ub(1)-A(lint,noe)*ub(2);
elseif cb=='dn'
f(lint)=f(lint)-A(lint,1)*ub(1);
elseif cb=='nd'
f(lint)=f(lint)-A(lint,noe)*ub(2);
end

A=A(lint,lint);
f=f(lint);
noeint=i2-i1+1;

% MEG on the linear system

un=A\f;

if cb=='dd'
un=[ub(1);un;ub(2)];
elseif cb=='dn'
un=[ub(1);un];
elseif cb=='nd'
un=[un;ub(2)];
end

return
