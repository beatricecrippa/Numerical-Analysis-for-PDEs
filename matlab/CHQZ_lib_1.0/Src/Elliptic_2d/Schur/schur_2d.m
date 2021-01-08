function [xy,un,param]=schur_2d(xa,xb,ya,yb,gam,...
          uex,uex_x,uex_y,ff,g,h,cb,nex,nx,ney,ny,gammax,gammay,param);
% SCHUR_2D   Numerical solution of the 2D b.v.p. -Delta u + gam u -Schur complement
%
%     -Delta u + gam u = f       x in Omega
%
%      + Dirichlet bc
%
% by  Schur Complement Matrix (Sect. 6.4.3, CHQZ3)
%  and by SEM Numerical Integration with LGL quadrature formulas.
%
%  [xy,un,param]=schur_2d(xa,xb,ya,yb,gam,...
%         uex,uex_x,uex_y,ff,nex,nx,ney,ny,gammax,gammay,param);
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
%        param(2)     : NOT USED
%        param(3)     : NOT USED
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
%        param(10)= tolerance for PCG stopping test
%        param(11)= maximum number of iterations for PCG stopping test
%        param(12:20)  NOT USED
%                   
% Output: xy = 2-indexes array wiht coordinates of 2D LGL mesh              
%         un = numerical solution
%         param(21) = (if param(3)==1)  iterations of PCG
%         param(22) = (if param(3)==1)  final residual 
%         param(23,24) = (if param(3)==2) extrema eigenvalues (max,min) of
%                  the preconditioned Schur complement
%         param(25,26,27) = (if param(4)==1 & param(3)==1) 
%                   errors on the exact solution:
%                   L^inf-norm, L2-norm, H1-norm
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

% Rigth Hand Side
f=ff(xy(:,1),xy(:,2)).*ww;
ub=ifro.*uex(xy(:,1),xy(:,2));

% Generation of local matrices and local r.h.s

[AGG,Amm,AGm,Lmm,LGG,Am,f]=schur_loc(ifro,nov,wx,dx,jacx,...
    wy,dy,jacy,nvli,gam,f,ub,param);

% Restrict r.h.s to internal nodes
f=f(lint);

noei=length(lint);
ifroi=ifro(lint);
xyi=[xy(lint,1),xy(lint,2)];
wwi=ww(lint);

% Generation of novg: its implements the 
%         restriction map R_{\Gamma_m} from \Gamma  to \Gamma_m (=
%         \partial\Omega_m\cap\Gamma)  (CHQZ, pag. 394)

ngamma=length(lgamma);
[novg]=cosnovg(xyi,noei,ifroi,lgamma,ldnov,novi,nvli);

%         Rgamma = matrix of size (ne, ngamma) (defined in CHQZ3, pag. 399)
%                  (R_\Gamma)_ij: =  1/n_j     if  x_j \in \Gamma_i
%                                 =  0         otherwise

[Rgamma]=cosrgam(novg,LGG,ne,ngamma);
fsigma=f(lgamma);

% compute f_sigma with formula (6.4.31) pag. 394 CHQZ3
for ie=1:ne
ngammam=length(LGG{ie});
agm=AGm{ie};amm=Amm{ie};lmm=Lmm{ie};lgg=LGG{ie};
floc=agm*(amm*f(novi(lmm,ie)));
fsigma(novg(lgg,ie))=fsigma(novg(lgg,ie))-floc;
end


u0=zeros(ngamma,1);
if (param(1)==1)
    
 % identity preconditioner 
 
[usigma,iter,res]=schur_pcg(u0, fsigma, param(10),param(11),param,...
    AGG,Amm,AGm,LGG,Am,nvli,novg);

elseif param(1)==2  % neumann-neumann preconditioner

% Computes the diagonal weighting matrix D relative to the interface
%  unknowns: 
%  D_ii=1/n_i   where n_i is the number of subdomains x_i belongs to.
%  D_m=D(novg(LGG_m,m))  is the restriction of D to Gamma_m 
%  (used in (6.4.40), pag. 397, CHQZ3 

D=partition(Rgamma);
[usigma,iter,res]=schur_pcg(u0, fsigma, param(10),param(11),param,...
    AGG,Amm,AGm,LGG,Am,nvli,novg,D);

elseif param(1)==3  % balancing neumann-neumann preconditioner

% Computes the diagonal weighting matrix D relative to the interface
%  unknowns: 
%  D_ii=1/n_i   where n_i is the number of subdomains x_i belongs to.
%  D_m=D(novg(LGG_m,m))  is the restriction of D to Gamma_m 
%  (used in (6.4.40), pag. 397, CHQZ3 

D=partition(Rgamma);

% Computes PSH= pseudoinverse of Sigma_H

[PSH]=pinv_sigma(AGG,Amm,AGm,LGG,novg,Rgamma);

[usigma,iter,res]=schur_pcg(u0, fsigma, param(10),param(11),param,...
    AGG,Amm,AGm,LGG,Am,nvli,novg,D,Rgamma,PSH);
end


param(21)=iter;
param(22)=res;
 
% computes local solutions in Omega_ie

[un]=local_solver(Amm,AGm,Lmm,LGG,novi,nvli,nov,novg,lint,lgamma,usigma,f,ub);


if (param(4)==1)
    % errors computing  on the exact solution
[err_inf,err_h1,err_l2]=errors_2d(x,wx,dx,xx,jacx,y,wy,dy,yy,jacy,...
xy,ww,nov,un,uex,uex_x,uex_y,param);
param(25:27)=[err_inf;err_h1;err_l2];
end

if (param(8)>0)
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
[hf,ha]=plot_sem(fig,command,nex,ney,x,xx,jacx,y,yy,jacy,xy,ww,nov,un,param);

end
