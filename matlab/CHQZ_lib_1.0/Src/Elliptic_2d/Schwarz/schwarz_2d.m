function [xy,un,param]=schwarz_2d(xa,xb,ya,yb,gam,...
          uex,uex_x,uex_y,ff,g,h,cb,nex,nx,ney,ny,gammax,gammay,param);
% SCHWARZ_2D   Numerical solution of the 2D boundary value problem
%
%     -Delta u + gam u = f       x in Omega
%
%      + Dirichlet bc
%
% by  SEM Numerical Integration with LGL quadrature formulas, with
%  additive Schwarz preconditioner with overlapping elements and
%   coarse mesh (Sect. 6.3.3, CHQZ3)
%
%  [xy,un,param]=schwarz_2d(xa,xb,ya,yb,gam,...
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
%        It works only if cb='dddd'
%        nex = number of elements (equally spaced) along x-direction
%        nx = polynomial degree in each element (the same in each element)
%               along x-direction
%        ney = number of elements (equally spaced) along y-direction
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
%        param(1) =1 P=I
%                 =2 P=Additive Schwarz with overlap and coarse mesh
%        param(2) = number of added layers for extending spectral elements
%                   inside additive Schwarz preconditioner.
%                   If param(2)==1, the preconditioner is P^{m}_{as,H} 
%                                    (minimum overlap), pag. 377 CHQZ3)
%                   If param(2)==2, the preconditioner is P^{s}_{as,H} 
%                                    (small overlap), pag. 377 CHQZ3)
%                   param(2) is a positive integer less than min(nx,ny)
%        param(3) = 1   PCG is used to solve linear system
%        param(3) = 2   PBiCGStab is used to solve linear system
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
%        param(10)= tolerance for PCG/PBCGSTAB stopping test
%        param(11)= maximum number of iterations for PCG/PBCGSTAB stopping test
%        param(12:20)  INTERNAL USE
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

% Matrix assembling

A=stiff_2d_se(npdx,nex,npdy,ney,nov,wx,dx,jacx,wy,dy,jacy);
M=spdiags(ww,0,noe,noe);

if gam ~=0
    A=A+gam*M;
end

% Rigth Hand Side
f=ff(xy(:,1),xy(:,2)).*ww;

% Dirichlet boundary conditions
ub=uex(xy(:,1),xy(:,2));

% 

Ab=A(lint,ldir);
A=A(lint,lint);
f=f(lint)-Ab*ub(ldir);
noei=length(lint);
wwi=ww(lint);

% extended elements construction

[nove,nvle]=cosnovenew(nx,nex,ny,ney,nov,ifro,param(2));

% Construction of local stiffness Q1 matrices on extended elements

[Aq1,Abq1,wwq1,linte,ldire,nove]=stiffq1(ifro,nov,xy,nove,nvle);

% Construction of Q1 stiffness matrix on the coarse grid.

[Ac,Acb,wwc,lista_coarse,noec,novc,lintc,ldirc]=stiffq1H(nx,nex,ny,ney,xy,nov,ifro);
r0t=matr0t(nx,ny,xy,nov,novc,noec,lista_coarse);

% Unity partition

nvl=[mn*ones(ne,1),npdx*ones(ne,1),npdy*ones(ne,1)];
[p_unity]=partition_e(nov,nvl,noe);

% PCG or PBCGstab with Schwarz preconditioner

u0=zeros(noei,1);
param(13:16)=[nx,ny,nex,ney];
if param(3)==1
[un,iter,res]=schwarz_pcg(u0, f, param, p_unity,...
    xy, ww, A, nov, noei, lint, x,wx, y,wy, xx,jacx,yy,jacy,...
    Aq1, wwq1,linte,nove,nvle,...
    Ac, Acb,wwc,r0t,lista_coarse,noec,novc,lintc,ldirc);
elseif param(3)==2
[un,iter,res]=schwarz_pbcgstab(u0, f, param, p_unity,...
    xy, ww, A, nov, noei, lint, x,wx, y,wy, xx,jacx,yy,jacy,...
    Aq1, wwq1,linte,nove,nvle,...
    Ac, Acb,wwc,r0t,lista_coarse,noec,novc,lintc,ldirc);
end
param(21:22)=[iter,res];
ub(lint)=un;
un=ub;

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
[ha]=plot_sem_2d(fig,command,nex,ney,x,xx,jacx,y,yy,jacy,xy,ww,nov,un,param(9));

end
