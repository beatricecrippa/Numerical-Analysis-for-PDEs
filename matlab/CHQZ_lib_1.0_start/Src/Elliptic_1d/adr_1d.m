function [xy,un,err_inf,err_l2,err_h1,der]=adr_1d(xa,xb,nu,beta,gam,...
          uex,uexx,ff,cb,ne,p,nx,param);
% ADR_1D   Numerical solution of the 1D boundary value problem (1.2.52)-(1.2.53), pag. 17, CHQZ2
%
%     -(nu u' + beta(x) u)'+gam u=f      xa < x < xb
%
%     + Dirichlet or Neumann bc
%
% by Galerkin Numerical Integration with LGL quadrature formulas.
%
%  [xy,un,err_inf,err_l2,err_h1,der]=adr_1d(xa,xb,nu,beta,gam,uex,uexx,...
%             ff,cb,ne,p,nx,param);
%
% Input: xa, xb = extrema of computational domain Omega=(xa,xb)
%      nu   = viscosity (constant>0)
%      beta = coefficient of first order term (beta=@(x)[beta(x)], 
%             with .*, .^, ./)
%      gam  = coefficient of zeroth order term (constant>0)
%      uex  = exact solution (uex=@(x)[uex(x)], with .*, .^, ./)
%      uexx = first derivative of exact solution (uexx=@(x)[uexx(x)], 
%             with .*, .^, ./)
%      ff  = r.h.s. solution (ff=@(x)[ff(x)], with .*, .^, ./)
%      cb = string containing boundary conditions, for example
%           cb='dn' imposes Dirichlet in xa and Neumann in xb 
%      ne = number of elements (equally spaced)
%      p = parameter : p=1 P_1 finite element, nx=1 
%                      p=2 P_2 finite element, nx=2 
%                      p=3 P_3 finite element, nx=3 
%                      p=4 spectral element, with polynomial degree nx, to be set 
%      nx = polynomial degree in each element (the same in each element)
%             to be set only if p=4, otherwise, nx=p;
%      param(1) = 1: compute errors (L^inf-norm, L2-norm, H1-norm)
%                      on the exact solution
%                 2: no  errors are computed
%      param(2) = 0: LG quadrature formulas with high precision degree are
%                      used to compute norms (exact norms)
%                 1: LGL quadrature formulas with npdx,npdy nodes are
%                      used to compute norms (discrete norms)
%                   (used only if param(1) == 1)
%      param(3) = number of nodes for high degree quadrature formula,
%                   (used only if param(2) == 0 & param(1) == 1)
%      param(4) = 0: absolute errors are computed
%                 1: relative errors are computed
%                   (used only if param(1) == 1)
%      param(5) = 0: not plot the solution
%                 1: plot the solution 
%      param(6) = number of nodes in each element for plotting interpolated solution
%                   (used only if param(5) ==1)
%
% Output: xy = 1D mesh
%         un = numerical solution 
%         err_inf = ||u_ex-u_N|| with respect to discrete maximum-norm
%         err_h1 = ||u_ex-u_N|| with respect to discrete H1-norm
%         err_l2 = ||u_ex-u_N|| with respect to discrete L2-norm
%         der    =  |(nu u'+ beta u)|_(x=xb) |
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

xy=0;un=0;err_inf=0;err_l2=0;err_h1=0;der=0;
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

b=zeros(noe,1)+beta(xy);

% Matrix assembling
A=adr_1d_se(npdx,ne,nov,nu,b,gam,wx,dx,jacx);

% Setting Dirichlet boundary conditions on the matrix

if cb(1)=='d'; A(1,:)=zeros(1,noe); A(1,1)=1; end
if cb(2)=='d'; A(noe,:)=zeros(1,noe); A(noe,noe)=1; end

% Right hand side
f=zeros(noe,1);
f=ff(xy).*ww;

% Setting Dirichlet boundary conditions on the matrix on r.h.s

if cb(1)=='d'; 
f(1)=uex(xy(1));
else
f(1)=f(1)+(-nu*uexx(xy(1))+b(1)*uex(xy(1))) ;
end
if cb(2)=='d'; 
f(noe)=uex(xy(noe));
else
f(noe)=f(noe)+nu*uexx(xy(noe))+b(noe)*uex(xy(noe));
end

% MEG on the linear system

un=A\f;


if param(1)==1

% Evaluates exact solution to compute errors

u=uex(xy);

% computes difference
err=u-un;
nq=param(3); fdq=param(2); err_type=param(4);

% ||u_ex-u_N|| with respect to discrete maximum-norm
if err_type==0
err_inf=norm(err,inf);
else
err_inf=norm(err,inf)/norm(u,inf);
end    

% ||u_ex-u_N|| with respect to H1 norm
% fdq=0 LG quadrature formula with nq nodes on each element
% fdq=1 LGL quadrature formula with npdx nodes on each element
[err_h1]=normah1_1d(fdq, nq, err_type, un, uex, uexx,...
 x, wx, dx, xx, jacx, xy,nov);

% ||u_ex-u_N|| with respect to L2 norm

[err_l2]=normal2_1d(fdq, nq, err_type, un, uex, ...
 x, wx, xx, jacx, xy,nov);



% compute the flux  der=(nu * u'+beta u)|(x=1)

ie=ne;
un1=dx(npdx,:)/jacx(ie)*un(nov(1:npdx,ie));

der=abs(nu*un1+b(noe)*un(noe)-(nu*uexx(xy(noe))+b(noe)*uex(xy(noe))));
end

if(param(5)==1)
    % interpolate numerical solution to give a "good" plot
    
    fig=figure(...,
    'Name','SEM-NI solution of -(nu u'' + beta(x) u)''+gam u=f     xa < x < xb',...
    'Visible','on');
    ha=plot_sem_1d(fig,ne,x,wx,xx,jacx,xy,nov,un,param(6));
end
return
