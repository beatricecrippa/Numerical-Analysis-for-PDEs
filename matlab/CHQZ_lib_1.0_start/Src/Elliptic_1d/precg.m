function [x,iter,res]=precg(A, b, x0, tol, itmax,P,preco)
% PRECG Preconditioned Conjugate Gradient method with FEM preconditioners
%
%  [x,iter,res]=precg(A, b, x0, tol, itmax,P,preco)
%
% preconditioned conjugate gradient
%
%    P^{-1}Ax=P^{-1}b with A and P SPD
%    The linear system P z=r is solved by Cholesky factorization. P is
%    factorized once.
%
% Input:  A: matrix (n x n)
%         b: right hand side (n x 1)
%         x0: intial datum (n x 1)
%         tol: tolerance for the stopping test
%         itmax: maximum number of iterations
%         preco : choice of the preconditioner:
%            1 : P=K_FE            (4.4.45), pag. 221, CHQZ2
%                A=K_GNI
%            2 : P=M_FE^{-1}K_FE   (4.4.46), pag. 221, CHQZ2
%                A=M_GNI^{-1}K_GNI
%            3 : P=M_FE,d^{-1}K_FE   (4.4.47), pag. 221, CHQZ2
%                A=M_GNI^{-1}K_GNI
%            4 : P=M_FE^{-1/2}K_FE M_FE^{-1/2}  (4.4.48), pag. 221, CHQZ2
%                A=M_GNI^{-1/2}K_GNI M_GNI^{-1/2}
%            5 : P=KM_FE<d^{-1/2}K_FE M_FE,d^{-1/2}_FE (4.4.49), pag. 221, CHQZ2
%                 A=M_GNI^{-1/2}K_GNI M_GNI^{-1/2}
%         P: preconditioner
%
% Output: x: solution (n x 1)
%         iter: iterations needed to satisfy stopping test
%         res: normalized 2-norm residual at the last iteration iter
%
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$


x=x0;
res=[];
r=b-A*x;
if(preco==0)
% 
z=r;
elseif(preco==1 | preco==4 | preco ==5)
% 
P1=chol(P);
z=P1\((P1)'\r);

end
rr0=z'*r;
normar0=sqrt(r'*r);
p=z;
err=normar0; iter=0;

while  err > tol & iter < itmax

v=A*p;
alpha=rr0/(p'*v);
x=x+alpha*p;
r=r-alpha*v;
if(preco==0)
z=r;
elseif(preco==1 | preco==4 | preco==5)
z=P1\((P1)'\r);
end
rr=z'*r;
beta=rr/rr0;
p=z+beta*p;
err=sqrt(r'*r)/normar0;
%res=[res;err];
iter=iter+1;
rr0=rr;
%fprintf('It. GC %d, err %13.6e \n',iter,err)
end
res=err;
if err>tol & iter>=itmax
fprintf('The method does not converge to prefixed tolerance in %d iterations\n',iter)
end
