function [xyz,un,D,param]=lap_3d(xa,xb,ya,yb,za,zb,gam,...
          uex,uex_x,uex_y,uex_z,ff,nex,nx,ney,ny,nez,nz,gammax,gammay,gammaz,...
          param);
% LAP_3D   Numerical solution of the 3D b.v.p. -Delta u + gam u = f
%
%     -Delta u + gam u = f       x in Omega
%
%      + Dirichlet bc
%
% by SEM Numerical Integration with LGL quadrature formulas.
%
%  [xy,un,D,param]=lap_3d(xa,xb,ya,yb,za,zb,gam,...
%         uex,uex_x,uex_y,uex_z,ff,nex,nx,ney,ny,nez,nz,gammax,gammay,gammaz,...
%         param);
%
% Input: xa= abscissa of either vertex V1 
%        xb= abscissa of either vertex V2 
%        ya= ordinate of either vertex V1 
%        yb= ordinate of either vertex V3 
%        za= ordinate of either vertex V1 
%        zb= ordinate of either vertex V5 
%        gam   = coefficient of zeroth order term (constant>=0)
%        uex  = exact solution (uex=@(x,y,z)[uex(x,y,z)], with .*, .^, ./)
%        uex_x  = exact first x-derivative (uex_x=@(x,y,z)[uex_x(x,y,z)], 
%                 with .*, .^, ./)
%        uex_y  = exact first y-derivative (uex_y=@(x,y,z)[uex_y(x,y,z)], 
%                 with .*, .^, ./)
%        uex_z  = exact first z-derivative (uex_z=@(x,y,z)[uex_z(x,y,z)], 
%                 with .*, .^, ./)
%        ff  = r.h.s. solution (ff=@(x,y,z)[ff(x,y,z)], with .*, .^, ./)
%        nex = number of elements (equally spaced) along x-direction
%        nx = polynomial degree in each element (the same in each element)
%               along x-direction
%        ney = number of elements (equally spaced) along y-direction
%        ny = polynomial degree in each element (the same in each element)
%               along y-direction
%        nez = number of elements (equally spaced) along z-direction
%        nz = polynomial degree in each element (the same in each element)
%               along z-direction
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
%        gammaz = column or row array of length nez-1. 
%               If the deomposition of Omega is not
%               uniform, gammay is the vector of position of interfaces between
%               spectral elements along z-direction. If gammay=[], uniform
%               decomposition is used.
%        param = array of parameters
%        param(1) = 1 : SEM-NI approach, see Sect. 5.3.3, CHQZ3
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
% Output: xyz = 3-indexes array wiht coordinates of 3D LGL mesh              
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
%         param(27,28) = (if param(3)==6) extrema generalized eigenvalues 
%                 (max,min) of
%                 A v = lambda M v, where A is the stiffness matrix and M is the
%                 mass matrix
%         param(29,30,31) = (if param(4)==1 & (param(3)==1 | param(3)==3)) 
%                   errors on the exact solution:
%                   L^inf-norm, L2-norm, H1-norm
%
%
%           Omega=(xa,xb) x (ya,yb) x (za,zb)
%
%             V8  _________ V7
%               /|        /|
%              / |       / |
%          V5 /________ /V6|
%             |  |      |  |
%             |V4|______|__| V3
%             |  /      | /               
%             | /       |/
%             |/________/
%           V1            V2
%                      
%
%
%       __________________________
%       |      |      |     |     |
%       |  3   |  6   |  9  | 12  |      Spectral elements
%       |      |      |     |     |      ordering in a plane yz then 
%       __________________________       planes at different x follow.
%       |      |      |     |     |
%       |  2   |  5   |  8  | 11  |
%       |      |      |     |     |
%       __________________________
%       |      |      |     |     |
%       |  1   |  4   |  7  | 10  |
%       |      |      |     |     |
%       __________________________
%

xyz=[]; un=[]; D=[]; 

npdx=nx+1; npdy=ny+1; npdz=nz+1;
mn=npdx*npdy; ldnov=mn*npdz; ne=nex*ney*nez;


[x,wx]=xwlgl(npdx); [dx]=derlgl(x,npdx);
[y,wy]=xwlgl(npdy); [dy]=derlgl(y,npdy);
[z,wz]=xwlgl(npdz); [dz]=derlgl(z,npdz);

% Mesh generation (with nov construction)

[xx,yy,zz,jacx,jacy,jacz,xyz,ww,ifro,nov]=mesh_3d(xa,xb,ya,yb,za,zb,nex,ney,nez,...
npdx,npdy,npdz,x,wx,y,wy,z,wz,gammax,gammay,gammaz);

% Matrix assembling 
noe=length(ww);
A=stiff_3d_se(npdx,nex,npdy,ney,npdz,nez,nov,wx,dx,jacx,wy,dy,jacy,wz,dz,jacz);
M=spdiags(ww,0,noe,noe);

if gam ~=0
    A=A+gam*M;
end

% Rigth Hand Side
f=ff(xyz(:,1),xyz(:,2),xyz(:,3)).*ww;



% Lists of internal, boundary, interface nodes. All lists are referred to
% global ordering on Omega
%
% lbor: list of boundary nodes
% lint: list of internal nodes 
% lintint: list of internal nodes which are internal to spectral elements 
% lgamma: list of internal nodes which are on the interface between spectral elements 



[ldir,lint,lintint,lgamma,ifro]=liste(ifro,nov);


% Setting Dirichlet boundary conditions on both matrix and r.h.s
ub=ifro.*uex(xyz(:,1),xyz(:,2),xyz(:,3));

f(lint)=f(lint)-A(lint,ldir)*ub(ldir);
A=A(lint,lint); M=M(lint,lint);
f=f(lint);
noei=length(lint);

%   selection on param(3)


if (param(3)==1)

    % Solve the linear system by a direct method

if param(2)==0
        % no reordering
        t_ord=0;
t=cputime; un=A\f;t_solv=cputime-t;

elseif param(2)==1
        % Cuthill-McKee ordering
t=cputime; pp=symrcm(A);t_ord=cputime-t;
t=cputime; un(pp)=A(pp,pp)\f(pp);t_solv=cputime-t;

elseif param(2)==2
    % Symmetric Minimum degree ordering
t=cputime; pp=symamd(A);t_ord=cputime-t;
t=cputime; un(pp)=A(pp,pp)\f(pp);t_solv=cputime-t;
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
[err_inf,err_h1,err_l2]=errors_3d(x,wx,dx,xx,jacx,y,wy,dy,yy,jacy,...
z,wz,dz,zz,jacz,...
xyz,ww,nov,un,uex,uex_x,uex_y,uex_z,param);
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
%[hf,ha]=plot_sem(fig,command,nex,ney,x,xx,jacx,y,yy,jacy,xy,ww,nov,un,param);

end
