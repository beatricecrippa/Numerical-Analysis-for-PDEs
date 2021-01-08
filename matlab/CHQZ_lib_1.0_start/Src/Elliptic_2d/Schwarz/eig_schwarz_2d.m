function [param]=eig_schwarz_2d(xa,xb,ya,yb,gam,...
          cb,nex,nx,ney,ny,gammax,gammay,param);
% EIG_SCHWARZ_2D   Eigenvalues computation for the matrix associated to 2D b.v.p. -Delta u + gam u = f + Dirichlet b.c. 
%
%
%     -Delta u + gam u = f       x in Omega
%
%      + Dirichlet bc
%
%  SEM Numerical Integration with LGL quadrature formulas, with
%  additive Schwarz preconditioner with overlapped elements and
%  coarse mesh (Sect. 6.3.3, CHQZ3)
%
%  [param]=eig_schwarz_2d(xa,xb,ya,yb,gam,...
%          nex,nx,ney,ny,gammax,gammay,param);
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
%        cb= a string of four elements, each for any side of \partial\Omega. 
%           cb(i) refers to side number "i" of \partial\Omega
%           'd' character stands for Dirichlet boundary conditions on that side
%           'n' character stands for Neumann boundary conditions on that side
%        It works only if cb='dddd'
%        nex = number of elements (equally spaced) along x-direction
%        nx = polynomial degree in each element (the same in each element)
%               along x-direction
%        ney = number of elements (equally spaced) along y-direction
%        ny = polynomial degree in each element (the same in each element)
%               along y-direction%        gammax = column or row array of length nex-1. 
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
%        param(1) =1 P=I
%                 =2 P=Additive Schwarz with overlap and coarse mesh
%        param(2) = number of added layers for extending spectral elements
%                   inside additive Schwarz preconditioner.
%                   If param(2)==1, the preconditioner is P^{m}_{as,H} 
%                                    (minimum overlap), pag. 377 CHQZ3)
%                   If param(2)==2, the preconditioner is P^{s}_{as,H} 
%                                    (small overlap), pag. 377 CHQZ3)
%                   param(2) is a positive integer less than min(nx,ny)
%        param(3)     : maximum number of iterations in eigs
%        param(4)     : tolerance for stopping test in eigs
%        param(5:20)  NOT USED or INTERNAL USE
%                   
% Output: 
%         param(11,12)  extrema eigenvalues (max,min) of
%                  the preconditioned stiffness matrix
%         param(13)  iterative condition number (|lam_max|/|lam_min|) of
%                  the preconditioned stiffness matrix
%
% References: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.
%             CHQZ3 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Evolution to Complex Geometries 
%                     and Applications to Fluid DynamicsSpectral Methods"
%                    Springer Verlag, Berlin Heidelberg New York, 2007.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

xy=[]; un=[];  

npdx=nx+1; npdy=ny+1; ldnov=npdx*npdy; mn=ldnov; ne=nex*ney;

[x,wx]=xwlgl(npdx); [dx]=derlgl(x,npdx);
[y,wy]=xwlgl(npdy); [dy]=derlgl(y,npdy);

% Generation of nov: its implements the 
%         restriction map R_m from \overline\Omega  to overline\Omega_m
%         (CHQZ, pag. 394)

[nov]=cosnov_2d(npdx,nex,npdy,ney);
noe=nov(ldnov,ne);

% Mesh generation

[xx,yy,jacx,jacy,xy,ww,ifro]=mesh_2d(xa,xb,ya,yb,cb,nex,ney,npdx,npdy,...
nov,x,wx,y,wy,gammax,gammay);


% Generation of lists of internal, boundary, interface nodes. All lists are referred to
% global ordering on Omega
%
% ldir: list of Dirichlet boundary nodes
% lint: list of internal nodes 
% lintint: list of internal nodes which are internal to spectral elements 
% lgamma: list of internal nodes which are on the interfaces between spectral elements 

[ldir,lint,lintint,lgamma,ifro]=liste(ifro,nov);

% Matrix assembling

A=stiff_2d_se(npdx,nex,npdy,ney,nov,wx,dx,jacx,wy,dy,jacy);
M=spdiags(ww,0,noe,noe);

if gam ~=0
    A=A+gam*M;
end


% 

A=A(lint,lint);
noei=length(lint);
wwi=ww(lint);

% extended elements construction

[nove,nvle]=cosnovenew(nx,nex,ny,ney,nov,ifro,param(2));

% Construction of local stiffness Q1 matrices on extended elements

[Aq1,Abq1,wwq1,linte,ldire,nove]=stiffq1(ifro,nov,xy,nove,nvle);

% Construction of stiffness Q1 matrix on the coarse grid.

[Ac,Acb,wwc,lista_coarse,noec,novc,lintc,ldirc]=stiffq1H(nx,nex,ny,ney,xy,nov,ifro);
r0t=matr0t(nx,ny,xy,nov,novc,noec,lista_coarse);

% Unity partition

nvl=[mn*ones(ne,1),npdx*ones(ne,1),npdy*ones(ne,1)];
[p_unity]=partition_e(nov,nvl,noe);
param(13:16)=[nx,ny,nex,ney];

% Matrix construction

Pas=zeros(noei,noei);
for i=1:noei
    r=A(:,i);
    [z]=precoasc(r,param,noei,lint,p_unity,...
    xy, ww, nov, x,wx,y,wy,xx,jacx,yy,jacy,...
    Aq1, wwq1, linte, nove, nvle,...
    Ac, Acb, wwc, r0t,lista_coarse,noec,novc,lintc,ldirc);
Pas(:,i)=z;
end

% Eigenvalues computation

if param(4)==0
param(4)=1000;
end
if param(5)==0
param(5)=1.d-12;
end
opts.issym=1;opts.disp=0;
opts.maxit=param(4); opts.tol=param(5);
lambda_max=eigs(Pas,1,'lm',opts);
lambda_min=eigs(Pas,1,'sm',opts);
param(11)=lambda_max;
param(12)=lambda_min;
param(13)=lambda_max/lambda_min;


return
