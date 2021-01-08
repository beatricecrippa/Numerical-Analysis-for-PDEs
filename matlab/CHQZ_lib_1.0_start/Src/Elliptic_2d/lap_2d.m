function [xy,un,D,param]=lap_2d(xa,xb,ya,yb,gam,...
          uex,uex_x,uex_y,ff,g,h,cb,nex,nx,ney,ny,gammax,gammay,param);
% LAP_2D   Numerical solution of the 2D b.v.p. -Delta u + gam u = f
%
%     -Delta u + gam u = f       x in Omega
%
%      + Dirichlet bc
%
% by SEM Numerical Integration with LGL quadrature formulas.
%
%  [xy,un,D,param]=lap_2d(xa,xb,ya,yb,gam,...
%         uex,uex_x,uex_y,ff,g,h,nex,nx,ney,ny,gammax,gammay,param);
%
%           Omega=(xa,xb) x (ya,yb)
%
%                side 3
%           V4 __________ V3
%             |          |
%             |          |
%     side 4  | Omega    |  side 2
%             |          |
%             |__________|
%           V1            V2
%                side 1
%
%
%       __________________________
%       |      |      |     |     |
%       |  3   |  6   |  9  | 12  |      Omega and spectral elements
%       |      |      |     |     |      ordering
%       __________________________
%       |      |      |     |     |
%       |  2   |  5   |  8  | 11  |
%       |      |      |     |     |
%       __________________________
%       |      |      |     |     |
%       |  1   |  4   |  7  | 10  |
%       |      |      |     |     |
%       __________________________
%
%
% Input: xa= abscissa of either vertex V1 or vertex V4
%        xb= abscissa of either vertex V2 or vertex V3
%        ya= ordinate of either vertex V1 or vertex V2
%        yb= ordinate of either vertex V3 or vertex V4
%        gam   = coefficient of zeroth order term (constant>=0)
%        uex  = exact solution (uex=@(x,y)[uex(x,y)], with .*, .^, ./)
%        uex_x  = exact first x-derivative (uex_x=@(x,y)[uex_x(x,y)], with .*, .^, ./)
%        uex_y  = exact first y-derivative (uex_y=@(x,y)[uex_y(x,y)], with .*, .^, ./)
%        ff  = r.h.s. solution (ff=@(x,y)[ff(x,y)], with .*, .^, ./)
%        g  = @(x,y)[....] function handle to the expression of Dirichlet
%              boundary data
%        h =@(x,y)[....] function handle to the expression of Neumann
%              boundary data. It is a vector of 4 functions, each for any
%              side.
%        cb= a string of four elements, each for any side of \partial\Omega. 
%           cb(i) refers to side number "i" of \partial\Omega
%           'd' character stands for Dirichlet boundary conditions on that side
%           'n' character stands for Neumann boundary conditions on that side
%        nex = number of elements along x-direction
%        nx = polynomial degree in each element (the same in each element)
%               along x-direction
%        ney = number of elements along y-direction
%        ny = polynomial degree in each element (the same in each element)
%               along y-direction
%        gammax = column or row array of length nex-1. 
%               If the deomposition of Omega is not
%               uniform, gammax is the vector of position of interfaces between
%               spectral elements along x-direction. If gammax=[], uniform
%               decomposition is used.
%        gammay = column or row array of length ney-1. 
%               If the deomposition of Omega is not
%               uniform, gammay is the vector of position of interfaces between
%               spectral elements along y-direction. If gammay=[], uniform
%               decomposition is used.
%        param = array of parameters
%        param(1) = 1 : SEM-NI approach, see Sect. 5.3.3, CHQZ3
%                 = 2 : Patching approach: See Sect. 5.13, CHQZ3
%        param(2) = 0 : no reordering of the matrix, see pag. 194, CHQZ2
%                   1 : Cuthill-McKee ordering (by symrcm function of
%                   Matlab)
%                   2 : Symmetric minimum degree ordering (by symamd function of
%                   Matlab)
%                   ONLY used if param(3)==1
%        param(3) = 1 : Solve the differential problem with a direct method 
%                   (Cholesky) and compute errors on the exact solution
%                   2 : Compute extrema eigenvalues of stiffness matrix
%                   3 : Compute interface Schur complement, solve the diff.
%                   problem by Schur complement
%                   4: Compute Schur complement and its extrema eigenvalues
%                   5: Compute and plot all the eigenvalues of A (stiffness
%                   matrix)
%                   6: Compute extrema eigenvalues of  M^{-1}A, i.e. 
%                      the generalized eigenvalues of  Av=lambda* M v,
%                      where M is the SEM-NI mass matrix
%        param(4) = 1: compute errors (L^inf-norm, L2-norm, H1-norm) 
%                      on the exact solution
%                   2: no  errors are computed
%        param(5) = 0: LG quadrature formulas with high precision degree are
%                      used to compute norms (exact norms)
%                   1: LGL quadrature formulas with npdx,npdy nodes are
%                      used to compute norms (discrete norms)
%                   (used only if param(4) == 1)
%        param(6) = number of nodes for high degree quadrature formula,
%                   (used only if param(5) == 0 & param(4) == 1)
%        param(7) = 0: absolute errors are computed
%                   1: relative errors are computed
%                   (used only if param(4) == 1)
%        param(8) = 0: not plot the solution
%                   1: plot the solution (mesh)
%                   2: plot the solution (surf)
%                   3: plot the solution (contour)
%        param(9) = number of nodes in each element and along each direction 
%                   for plotting the solution
%                   (used only if param(8) > 0)
%                   
% Output: xy = 2-indexes array wiht coordinates of 2D LGL mesh              
%         un = numerical solution
%         D  = (if param(3)==5)  eigenvalues of A
%         param(21) = (if param(3)==1)  cputime for reordering stiffness
%                  matrix befor solving
%         param(22) = (if param(3)==1)  cputime for solving linear system 
%                  by Choleski factorization
%         param(23,24) = (if param(3)==2) extrema eigenvalues (max,min) of
%                  stiffness matrix A
%         param(25,26) = (if param(3)==4) extrema eigenvalues (max,min) of
%                  interface Schur complement
%         param(27,28) = (if param(3)==6) extrema generalized eigenvalues (max,min) of
%                   A v = lambda M v, where A is the stiffness matrix and M is the
%                   mass matrix
%         param(29,30,31) = (if param(4)==1 & (param(3)==1 | param(3)==3)) 
%                   errors on the exact solution:
%                   L^inf-norm, L2-norm, H1-norm
%
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$



xy=[]; un=[]; D=[]; 

npdx=nx+1; npdy=ny+1; ldnov=npdx*npdy; mn=ldnov; ne=nex*ney;

[x,wx]=xwlgl(npdx); [dx]=derlgl(x,npdx);
[y,wy]=xwlgl(npdy); [dy]=derlgl(y,npdy);
%
% nov construction
%
[nov]=cosnov_2d(npdx,nex,npdy,ney);
noe=nov(ldnov,ne);

% Mesh generation

[xx,yy,jacx,jacy,xy,ww,ifro]=mesh_2d(xa,xb,ya,yb,cb,nex,ney,npdx,npdy,...
nov,x,wx,y,wy,gammax,gammay);

% Matrix assembling 

A=stiff_2d_se(npdx,nex,npdy,ney,nov,wx,dx,jacx,wy,dy,jacy);
M=spdiags(ww,0,noe,noe);

if gam ~=0
    A=A+gam*M;
end

% Rigth Hand Side
f=ff(xy(:,1),xy(:,2)).*ww;

% Neumann boundary condition
[f]=neumannbc(f,h,jacx,jacy,wx,wy,xy,ifro,nov);

%
% Patching approach 
% 
if (param(1)==2)
[A,f]=patch_se(A,f,ifro,nov,dx,jacx,dy,jacy);
end

% Lists of internal, boundary, interface nodes. All lists are referred to
% global ordering on Omega
%
% lbor: list of Dirichlet boundary nodes
% lint: list of internal nodes 
% lintint: list of internal nodes which are internal to spectral elements 
% lgamma: list of internal nodes which are on the interface between spectral elements 



[ldir,lint,lintint,lgamma,ifro]=liste(ifro,nov);


% Setting Dirichlet boundary conditions on both matrix and r.h.s
ub=zeros(noe,1);
for i=1:noe
if(ifro(i)==1)
ub(i)=uex(xy(i,1),xy(i,2));
end
end

if length(ldir)>0
f(lint)=f(lint)-A(lint,ldir)*ub(ldir);
A=A(lint,lint); M=M(lint,lint);
f=f(lint);
noei=length(lint);
else
    noei=noe;
end

%   selection on param(3)


if (param(3)==1)

if (param(1)==1)
    % SEM-NI discretization A is symmetric
    % Solve the linear system by Cholesky factorization

if param(2)==0
        % no reordering
        t_ord=0;
t=cputime; L=chol(A); un=L\(L'\f);t_solv=cputime-t;

elseif param(2)==1
        % Cuthill-McKee ordering
t=cputime; pp=symrcm(A);t_ord=cputime-t;
t=cputime; L=chol(A(pp,pp)); un(pp)=L\(L'\f(pp));t_solv=cputime-t;

elseif param(2)==2
    % Symmetric Minimum degree ordering
t=cputime; pp=symamd(A);t_ord=cputime-t;
t=cputime; L=chol(A(pp,pp)); un(pp)=L\(L'\f(pp));t_solv=cputime-t;
end
elseif (param(1)==2)
    % Patching discretization A is not symmetric
    % Solve the linear system by LU factorization

if param(2)==0
        % no reordering
        t_ord=0;
t=cputime; un=A\f;t_solv=cputime-t;

elseif param(2)==1
        % Cuthill-McKee ordering
t=cputime; pp=symrcm(A);t_ord=cputime-t;
t=cputime; [L,U,PP]=lu(A(pp,pp)); un(pp)=U\(L\(PP*f(pp)));t_solv=cputime-t;

elseif param(2)==2
    % Symmetric Minimum degree ordering
t=cputime; pp=symamd(A);t_ord=cputime-t;
t=cputime; [L,U,PP]=lu(A(pp,pp)); un(pp)=U\(L\(PP*f(pp)));t_solv=cputime-t;
end
end

param(21)=t_ord;
param(22)=t_solv;
 
%
ub(lint)=un;
un=ub;



elseif(param(3)==2)
    
% Compute extrema eigenvalues of the stiffness matrix
keigs=1;
if(param(1)==1)
    % Symmetric problem
opts.issym=1;opts.disp=0;opts.maxit=1000; opts.tol=1.d-12;
lambda_max=eigs(A,keigs,'lm',opts);
lambda_min=eigs(A,1,'sm',opts);
else
opts.issym=0;opts.disp=0;opts.maxit=1000; opts.tol=1.d-12;
lambda_max=eigs(A,keigs,'lm',opts);
lambda_min=eigs(A,1,'sm',opts);
end
param(23)=lambda_max;
param(24)=lambda_min;


elseif(param(3)==3)

% Compute interface Schur complement, solve the differential problem by 
% the Schur complement

Sigma=A(lgamma,lgamma)-A(lgamma,lintint)*inv(A(lintint,lintint))*A(lintint,lgamma);

% Solve  Schur linear system and internal linear systems by \ of matlab

fsigma=f(lgamma)-A(lgamma,lintint)*inv(A(lintint,lintint))*f(lintint);
xsigma=Sigma\fsigma;
xint=A(lintint,lintint)\(f(lintint)-A(lintint,lgamma)*xsigma);
u2=zeros(noei,1);
u2(lgamma)=xsigma;
u2(lintint)=xint;
ub(lint)=u2;
un=ub;


elseif(param(3)==4)

% Compute interface Schur complement and its extrema eigenvalues 

Sigma=A(lgamma,lgamma)-A(lgamma,lintint)*inv(A(lintint,lintint))*A(lintint,lgamma);

clear A; clear f;
opts.issym=1;opts.disp=0;opts.maxit=1000; opts.tol=1.d-12;
lambda_max=eigs(Sigma,1,'lm',opts);
lambda_min=eigs(Sigma,1,'sm',opts);
param(25)=lambda_max;
param(26)=lambda_min;


elseif(param(3)==5)
% compute and plot eigenvalues of A
D=eig(full(A));
fig=figure(...,
    'Name','Spectrum of SEM-NI stiffness matrix',...
    'Visible','on')
plot(real(D),imag(D),'o')
set(gca,'Fontname','Times','Fontsize',18)
xlabel('Re'); ylabel('Im')
if param(1)==1
set(gca,'Fontname','Times','Fontsize',18)
title(['N=',num2str(nx),' H=',num2str(2/nex), ' SEM-NI'])
else
set(gca,'Fontname','Times','Fontsize',18)
title(['N=',num2str(nx),' H=',num2str(2/nex), ' Patching'])
end


elseif(param(3)==6)
    keigs=1;
% compute extrema eigenvalues of  M^{-1}A, i.e. 
%  the generalized eigenvalues of  Av=lambda* M v
if(param(1)==1)
opts.issym=1;opts.disp=0;opts.maxit=1000; opts.tol=1.d-12;
lambda_max=eigs(A,M,keigs,'lm',opts);
lambda_min=eigs(A,M,keigs,'sm',opts);
param(27)=lambda_max;
param(28)=lambda_min;

elseif(param(2)==2)
    disp('No eigenvalues are computed')
end

end

if (param(4)==1 & (param(3)==1 | param(3)==3))
    % compute errors on the exact solution
[err_inf,err_h1,err_l2]=errors_2d(x,wx,dx,xx,jacx,y,wy,dy,yy,jacy,...
xy,ww,nov,un,uex,uex_x,uex_y,param);
param(29:31)=[err_inf;err_h1;err_l2];
end

if (param(8)>0 & (param(3)==1 | param(3)==3))
    if param(8)==1
        command='mesh';
    elseif param(8)==2
        command='surf';
    elseif param(8)==3;
        command='contour'
    end
% plot the solution
fig=figure(...,
    'Name','SEM-NI solution of -Delta u + gam u=f, Dirichlet bc',...
    'Visible','on');
[ha]=plot_sem_2d(fig,command,nex,ney,x,xx,jacx,y,yy,jacy,xy,ww,nov,un,param(9));

end

return
