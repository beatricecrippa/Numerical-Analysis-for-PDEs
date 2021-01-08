function [x,iter,res]=schur_pcg(x0, b, tol, maxit,param,...
    AGG,Amm,AGm,LGG,Am,nvli,novg,D,Rgamma,PSH)
% SCHUR_PCG   Preconditioned conjugate gradient to solve the Schur complement system
%      Preconditioned conjugate gradient to solve the Schur complement system
%
%  [x,iter,res]=schur_pcg(x0, b, tol,maxit,param,...
%               AGG,Amm,AGm,LGG,Am,nvl,novg,D,Rgamma,PSH)
%
% Input: x0 = column array for intial guess
%        b = r.h.s for PCG
%        tol = tolerance for PCG stopping test
%        maxit = maximum number of iterations for PCG stopping test
%        param = array of parameters
%        AGG, Amm, AGm, LGG, Am = cells produced in schur_local.m.
%              They contain local matrices and lists 
%        nvli = column array. nvli(ie) is the number of nodes of \Omega_ie
%        internal to Omega.
%         novg(i,ie)= 0 if node x_i of Omega_ie does not belong to Gamma
%                    = j if node x_i of Omega_ie is the node j of Gamma 
%        D = diagonal weighting matrix relative to the interface
%            unknowns, set in partition.m
%        Rgamma = matrix of size (ne, ngamma) (defined in CHQZ3, pag. 399)
%                  set in cosrgam.m
%        PSH = pseudoinverse of Sigma_H, set in pinv_sigma.m
%
% Output: x = column array of PCG solution
%         iter = number of iterations required to satisfy stopping test
%                ||r^(k)||/||b|| < tol
%         res =  value of ||r^(k)||/||b|| at last iteration.
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

ne=length(nvli);
x=x0;
bb=norm(b);
res=[];
v=schur_mxv(x,AGG,Amm,AGm,LGG,novg,ne);
r=b-v;
if(param(1)==1)
z=r;
elseif(param(1)==2)
z=schur_preconnl(r,ne,nvli,novg,D,LGG,Am);
elseif(param(1)==3)
z=schur_precobnn(r,ne,nvli,novg,D,LGG,Am,Rgamma,PSH,AGG,Amm,AGm);
end
rr0=z'*r;
normar0=sqrt(r'*r);
p=z;
err=normar0/bb; iter=0;


while  err > tol & iter < maxit

%v=A*p;
v=schur_mxv(p,AGG,Amm,AGm,LGG,novg,ne);
alpha=(r'*z)/(p'*v);
x=x+alpha*p;
r=r-alpha*v;
if(param(1)==1)
z=r;
elseif(param(1)==2)
z=schur_preconnl(r,ne,nvli,novg,D,LGG,Am);
elseif(param(1)==3)
z=schur_precobnn(r,ne,nvli,novg,D,LGG,Am,Rgamma,PSH,AGG,Amm,AGm);
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
if err>tol & iter>=maxit
fprintf('No convergence is reached in  %d iterations \n',iter)
end

return
