function [u,iter,err]=pbcgstab_vdv(A, b, u, tol, maxit ,P,preco);
% PBCGSTAB_VDV Preconditioned BiCGStab method with FEM preconditioners
%
%  [x,iter,res]=pbcgstab_vdv(A, b, x0, tol, itmax,P,preco)
%
% preconditioned BiCGStab method
%
%    P^{-1}Ax=P^{-1}b with A and P SPD
%    The linear system P z=r is solved by LU factorization. P is
%    factorized once.
%
% Input:  A: matrix (n x n)
%         b: right hand side (n x 1)
%         x0: intial datum (n x 1)
%         tol: tolerance for the stopping test
%         itmax: maximum number of iterations
%         preco : choice of the preconditioner:
%           1 : P=K_FE            (4.4.45), pag. 221, CHQZ2
%               A=K_GNI
%           2 : P=M_FE^{-1}K_FE   (4.4.46), pag. 221, CHQZ2
%               A=M_GNI^{-1}K_GNI
%           3 : P=M_FE,d^{-1}K_FE   (4.4.47), pag. 221, CHQZ2
%               A=M_GNI^{-1}K_GNI
%           4 : P=M_FE^{-1/2}K_FE M_FE^{-1/2}  (4.4.48), pag. 221, CHQZ2
%               A=M_GNI^{-1/2}K_GNI M_GNI^{-1/2}
%           5 : P=KM_FE<d^{-1/2}K_FE M_FE,d^{-1/2}_FE (4.4.49), pag. 221, CHQZ2
%               A=M_GNI^{-1/2}K_GNI M_GNI^{-1/2}
%         P: preconditioner
%
% Output: x: solution (n x 1)
%         iter: iterations needed to satisfy stopping test
%         res: normalized 2-norm residual at the last iteration iter

%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

r=b-A*u;
rtilde0=r;
normar0=r'*r;
iter=0; err=1;
alpha=1;rho=1;rhom1=1;w=1;
 if preco~=0
 [L,U,PP]=lu(P);
 end

while iter< maxit & err> tol
rho=rtilde0'*r;
if iter==0
p=r;
else
beta=(rho/rhom1)*(alpha/w);
p=r+beta*(p-w*v);
end
if preco==0
pc=p;
else
pc=U\(L\(PP*p));
end

v=A*pc;
alpha=rho/(v'*rtilde0);
s=r-alpha*v;
if ( norm(s) < tol ),                          % early convergence check
        u = u + alpha*pc;
        err = norm( s ) / normar0;
        break;
     end
if preco==0
sc=s;
else
sc=U\(L\(PP*s));
end
t=A*sc;
w=(t'*s)/(t'*t);
u=u+alpha*pc+w*sc;
r=s-w*t;
rhom1=rho;
iter=iter+1;
err=sqrt(r'*r)/normar0;
end
