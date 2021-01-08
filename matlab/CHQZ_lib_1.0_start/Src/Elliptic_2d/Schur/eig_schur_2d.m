function [param]=eig_schur_2d(xa,xb,ya,yb,gam,...
          cb,nex,nx,ney,ny,gammax,gammay,param);
% EIG_SCHUR_2D   Eigenvalues computation Schur complement matrix
%
% Eigenvalues computation for the matrix associated to 
%                the 2D boundary value problem
%
%     -Delta u + gam u = f       x in Omega
%
%      + Dirichlet bc
%
% by  Schur Complement Matrix (Sect. 6.4.3, CHQZ3)
%  and by SEM Numerical Integration with LGL quadrature formulas.
%
%  [param]=eig_schur_2d(xa,xb,ya,yb,gam,...
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
%           It works only if cb='dddd';
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
%        param(1) = 1 : No preconditioner (P=I)
%                   2 : P= Neumann-Neumann (Alg. 6.4.3, CHQZ3)
%                   3 : P= balancing Neumann-Neumann (Alg. 6.4.4, CHQZ3)
%        param(2)     : maximum number of iterations in eigs
%        param(3)     : tolerance for stopping test in eigs
%                   
% Output: 
%         param(11,12)  extrema eigenvalues (max,min) of
%                  the preconditioned Schur complement
%         param(13)  iterative condition number (|lam_max|/|lam_min|) of 
%                  the preconditioned Schur complement
%
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
[novi,nvli]=cosnovi(nov,ifro,lint);


noei=length(lint);
ifroi=ifro(lint);
xyi=[xy(lint,1),xy(lint,2)];
wwi=ww(lint);

% Generation of novg: its implements the 
%         restriction map R_{\Gamma_m} from \Gamma  to \Gamma_m (=
%         \partial\Omega_m\cap\Gamma)  (CHQZ, pag. 394)

ngamma=length(lgamma);
[novg]=cosnovg(xyi,noei,ifroi,lgamma,ldnov,novi,nvli);

% Construction of LGG, to set Rgamma

LGG=cell(ne,1);
for ie=1:ne
ifro_l=ifro(nov(1:mn,ie));
[lbor,lint]=liste1(ifro_l);
ifroi=ifro_l(lint);
[lbor,lint,lintint,lg]=liste1(ifroi);
LGG(ie)={lg};
end

%         Rgamma = matrix of size (ne, ngamma) (defined in CHQZ3, pag. 399)
%                  (R_\Gamma)_ij: =  1/n_j     if  x_j \in \Gamma_i
%                                 =  0         otherwise

[Rgamma]=cosrgam(novg,LGG,ne,ngamma);

% Computes the diagonal weighting matrix D relative to the interface
%  unknowns: 
%  D_ii=1/n_i   where n_i is the number of subdomains x_i belongs to.
%  D_m=D(novg(LGG_m,m))  is the restriction of D to Gamma_m 
%  (used in (6.4.40), pag. 397, CHQZ3 


if (param(1)==1)
 D=[];
else
 D=partition(Rgamma);
end

% Generation of both Schur complement matrix and preconditioner

[Sigma,PNN]=schur_matrix(ifro,nov,wx,dx,jacx,...
    wy,dy,jacy,nvli,gam,novg,lint,lgamma,D,Rgamma,param);


if(param(1)~=1) % Neumann Neumann o balancing Neumann-Neumann
    Sigma=PNN*Sigma;
end

if param(2)==0 
param(2)=1000;
end
if param(3)==0
param(3)=1.d-12;
end
opts.issym=1;opts.disp=0;
opts.maxit=param(2); opts.tol=param(3);
lambda_max=eigs(Sigma,1,'lm',opts);
lambda_min=eigs(Sigma,1,'sm',opts);
param(11)=lambda_max;
param(12)=lambda_min;
param(13)=lambda_max/lambda_min;


return
