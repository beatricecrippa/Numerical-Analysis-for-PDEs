function [xy,un,A,M,AFE,MFE,MFEd,d,kappa,param]=ellprecofem_1d(xa,xb,nu,...
    beta, gam,ff,cb,ub,ne,nx,param);
% ELLPRECOFEM_1D  FEM-preconditioned SEM-appx of 1D second order b.v.p.
%
%   Numerical solution (or eigenvalues compuation) 
%          of the 1D boundary value problem
%
%     -nu u''+ beta*u' + gam u =f      xa < x < xb
%      Dirichlet or Neumann boundary conditions
%    
% by Galerkin Numerical Integration with LGL quadrature formulas
%    preconditioner by FEM matrices
%
% Useful to reproduce results on iterative condition numbers (1D case) 
% published in Ch4. of CHQZ2.
%
%    
%  [xy,un,A,M,AFE,MFE,MFEd,d,kappa,param]=ellprecofem_1d(xa,xb,nu,...
%        beta, gam,ff,cb,ub,ne,nx,param);
%
% Input: xa, xb = extrema of computational domain Omega=(xa,xb)
%      nu   = viscosity (constant>0)
%      beta  = coefficient of first order term (constant)
%      gam   = coefficient of zeroth order term (constant>=0)
%      ff  = r.h.s. solution (ff=@(x)[ff(x)], with .*, .^, ./)
%      cb = string containing boundary conditions, for example
%           cb='dn' imposes Dirichlet in xa and Neumann in xb,
%           possible b.c.: Dirichlet, Neumann
%      ub = 2 components array with Dirichlet, Neumann (pure derivative) 
%           data in xa and in xb
%      ne = number of elements (equally spaced)
%      nx = polynomial degree in each element (the same in each element)
%             to be set only if p=4, otherwise, nx=p;
%
%      param = parameters array :
%    param(1)= choice of the preconditioner:
%             0 : P=I
%                 A=K_GNI
%             1 : P=K_FE            (4.4.45), pag. 221, CHQZ2
%                 A=K_GNI
%             2 : P=M_FE^{-1}K_FE   (4.4.46), pag. 221, CHQZ2
%                 A=M_GNI^{-1}K_GNI
%             3 : P=M_FE,d^{-1}K_FE   (4.4.47), pag. 221, CHQZ2
%                 A=M_GNI^{-1}K_GNI
%             4 : P=M_FE^{-1/2}K_FE M_FE^{-1/2}  (4.4.48), pag. 221, CHQZ2
%                 A=M_GNI^{-1/2}K_GNI M_GNI^{-1/2}
%             5 : P=KM_FE<d^{-1/2}K_FE M_FE,d^{-1/2}_FE (4.4.49), pag. 221, CHQZ2
%                 A=M_GNI^{-1/2}K_GNI M_GNI^{-1/2}
%    param(2)= choice of the solver:
%             1 : Conjugate Gradient, only for symmetric preconditioners:
%                 usable only if param(1)=1, 4, 5
%             2 : Bi-CGStab (VDV), usable for any choice of param(1)
%            11 : eigenvalues computation, to have iterative condition 
%                 number 
%                 kapit(P^{-1}A)=|lambda_max(P^{-1}A)|/|lambda_min(P^{-1}A)| 
%            12 : singular values computation, to have 2-norm condition 
%                 number
%                 k_2(P^{-1}A)=sigma_max(P^{-1}A)/sigma_min(P^{-1}A)
%    param(3) = tolerance for the solver: 
%               stopping test is: ||r^k||_2/||r^0||_2 < tol
%               (not used if param(2)=11 | param(2)=12)
%    param(4) = maximum iteration number for the solver
%               (not used if param(2)=11 | param(2)=12)
%
%
% Output: xy = 1D mesh
%         un = numerical solution 
%         A = Stiffness matrix K_GNI (formula (4.4.43), pag. 219 CHQZ2)
%         M = Mass matrix K_GNI (formula (4.4.43), pag. 219 CHQZ2)
%         AFE = FEM Stiffness matrix K_FE (Tab. 4.6, pag. 221, CHQZ2)
%         MFE = FEM Mass matrix (Tab. 4.6, pag. 221, CHQZ2)
%         MFE = FEM Mass matrix with numerical integration 
%                      (Tab. 4.6, pag. 221, CHQZ2)
%         d = eigenvalues or singular values of the preconditioned matrix
%         kappa = iterative or 2-norm condition number (if param(2)=11, 12)
%         param(11) = iterations required by the solver
%         param(12) = residual at the last iterations
%         
%
% References: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.
%             VDV = van der Vorst, Henk A.,
%                   "Iterative Krylov methods for large linear systems"
%                   Cambridge Monographs on Applied and Computational Mathematics, 13. 
%                   Cambridge University Press, Cambridge,  2003.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

% check if the choice for the  preconditioner and the solver are  compatible
xy=[];un=[];A=[];M=[];AFE=[];MFE=[];MFEd=[];d=[];kappa=1;

if param(2)==1
    if(param(1)==2 | param(1)==3)
fprintf(' The choice for the  preconditioner and the solver are not compatible\n');
return
    end
end


npdx=nx+1;

[x,wx]=xwlgl(npdx); dx=derlgl(x,npdx);

% nov
ldnov=npdx; nov=zeros(ldnov,ne);
[nov]=cosnov_1d(npdx,ne,nov);
noe=nov(npdx,ne);

% Uniform decomposition of  Omega in ne elements

[xx,jacx,xy,ww]=mesh_1d(xa,xb,ne,npdx,nov,x,wx);

% ww contains the diagonal of the mass matrix

% G-NI  Stiffness and Mass Matrices assembling
if beta==0
A=stiff_1d_se(npdx,ne,nov,wx,dx,jacx); 
A=nu*A;
else
A=ad_1d_se(npdx,ne,nov,nu,beta,wx,dx,jacx); 
end    
M=spdiags(ww,0,noe,noe);
if gam ~=0
    A=A+gam*M;
end
% Right hand side
f=(ff(xy)).*ww;

% Setting Dirichlet boundary conditions on the matrix on r.h.s
if cb(1)=='d'
f(1)=ub(1); i1=2;
elseif cb(1)=='n'
f(1)=f(1)-nu*ub(1);i1=1;
end
if cb(2)=='d'; 
f(noe)=ub(2);i2=noe-1;
elseif cb(2)=='n'
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

f=f(lint);
noeint=i2-i1+1;

if param(1)~=0
% FEM  Stiffness and Mass Matrices assembling
[AFE,MFE,MFEd]=femp1_preco_1d(npdx,ne,nov,jacx,xy,nu,beta,gam,param);

% Setting Dirichlet boundary conditions on the preconditioned stiffness matrix

AFE=AFE(lint,lint);
MFE=MFE(lint,lint);
MFEd=MFEd(lint,lint);
end
A=A(lint,lint);
M=M(lint,lint);

% Pre solver
if param(1)==0
    P=speye(noeint);
    A1=A;f1=f;
elseif param(1)==1
    P=AFE;
    A1=A;f1=f;
elseif param(1)==2 
    P=M*inv(MFE)*AFE;
    A1=A;f1=f;
elseif param(1)==3
    P=M*inv(MFEd)*AFE;
    A1=A;f1=f;
elseif param(1)==4 
    MFE1=full(MFE);MFE1=sqrtm(MFE1);MFE1=inv(MFE1);
    M1=inv(sqrt(M));
    P=MFE1'*AFE*MFE1;
    A1=M1*A*M1;
    f1=M1*f;
elseif param(1)==5   
    MFEd1=inv(sqrt(MFEd)); 
    M1=inv(sqrt(M));
    P=MFEd1*AFE*MFEd1;
    A1=M1*A*M1;
    f1=M1*f;
end

% solve the linear system

u0=zeros(noeint,1);

if param(2)==1
    % CG to solve linear system
[un,iter,res]=precg(A1, f1, u0, param(3), param(4),P,param(1));
if param(1)==4 | param(1)==5
    un=M1*un;
end
if cb=='dd'
un=[ub(1);un;ub(2)];
elseif cb=='dn'
un=[ub(1);un];
elseif cb=='nd'
un=[un;ub(2)];
end
param(11:12)=[iter,res];

elseif param(2)==2
    % Bicgstab to solve linear system
[un,iter,res]=pbcgstab_vdv(A1, f1, u0, param(3), param(4),P,param(1));

if param(1)==4 | param(1)==5
    un=M1*un;
end
if cb=='dd'
un=[ub(1);un;ub(2)];
elseif cb=='dn'
un=[ub(1);un];
elseif cb=='nd'
un=[un;ub(2)];
end
param(11:12)=[iter,res];

elseif param(2)==11
% compute eigenvalues and the iterative condition number
C=inv(P)*A1;
OPTS.disp=0;
d1=eigs(C,3,'LM',OPTS);d2=eigs(C,3,'SM',OPTS);
kappa=max(abs(d1))/min(abs(d2));
elseif param(2)==12
% compute singular values and the 2-norm condition number
C=inv(P)*A1;
OPTS.disp=0;
d1=svds(C,3,'L',OPTS);d2=svds(C,3,0,OPTS);
kappa=max(d1)/min(d2);
end
return    
    
    
