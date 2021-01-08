function [x,u,err]=stat_scalar_hyp(xa,xb,beta,f,uex,ul,nx,param)
% STAT_SCALAR_HYP Numerical solution of a stationary scalar linear hyperbolic equation,
%     formula (3.7.1a), pag. 145, CHQZ2
%
%      beta du/dx=f(x)         in (xa,xb) 
%      inflow cond.: u(xa)=ul
%
%  Space discretization: 1. Galerkin-Numerical Integration with LGL formulas,
%                        2. Legendre collocation (strong form),
%                        3. penalty,
%
%
%  [x,u,err]=stat_scalar_hyp(xa,xb,beta,f,uex,ul,nx,param);
%
% Input: xa,xb  = extrema of space domain Omega=(xa,xb)
%        beta   = transport coefficient
%        uex    = exact solution (uex=@(x)[uex(x)], with .*, .^, ./)
%        ul    = inflow condition (scalar datum)
%        nx     = spectral polynomial degree, the number of LGL points is
%                 nx+1
%        param(1)= space disc. scheme: 1=strong (3.7.2), pag. 146 CHQZ2, 
%                                     2=GNI+ weak b.c. (3.7.5), pag. 147
%                                     3=penalty (3.7.2a)+(3.7.9), pag. 148
%        param(2)=tau : penalty parameter (used only if param(1)==3)
%
% Output: u = solution 
%         err= L_\infty norm error between numerical and exact solution 
%
% Possible input data:
% xa=-1;xb=1;
% beta=1.5;
% param(1)=2;param(2)=0.75
% nx=8;
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$



npdx=nx+1;
[x,wx]=xwlgl(npdx); 
[dx]=derlgl(x,npdx);
jac=(xb-xa)*5.d-1;
x=x*jac+(xa+xb)*5.d-1;


tau=param(2)/(jac*wx(1));

% Construction of matrix A

if param(1)==1
% strong collocation (3.7.2)
A=beta*diag(wx)*dx;
A(1,:)=zeros(1,npdx); A(1,1)=1;


elseif param(1)==2
% GNI+weak b.c. (3.7.5)
A=-beta*dx'*diag(wx);
A(npdx,npdx)=A(npdx,npdx)+beta;

elseif param(1)==3
% penalty (3.7.2a) +(3.7.9)
A=beta*diag(wx)*dx;
A(1,1)=A(1,1)+tau*beta*wx(1)*jac;

end

[L,U,P]=lu(A);

%right hand side
b=f(x).*wx*jac;


% inflow condition
if param(1)==1
b(1)=ul;
elseif param(1)==2
b(1)=b(1)+beta*ul;
elseif param(1)==3
b(1)=b(1)+tau*beta*ul*wx(1)*jac;
end

% linear system solution
u=U\(L\(P*b));


% error  

ue=uex(x);
err=norm(u-ue);

return
