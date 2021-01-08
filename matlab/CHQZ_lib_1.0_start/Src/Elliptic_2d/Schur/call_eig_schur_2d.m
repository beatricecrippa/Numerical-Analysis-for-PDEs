% CALL_EIG_SCHUR_2D Script for pre and post processing eig_schur_2d.m
%
%


%   Written by Paola Gervasio
%   $Date: 2007/04/01$

xa=-1;xb=1;   % Omega=(xa,xb) x (ya,yb)
ya=-1;yb=1;
gam=0;
cb='dddd';     %  eig_schur_2d works only if cb='dddd';
    param=zeros(20,1);  
    param(1)=3;    % 1:P=I, 2: P=NN, 3: P=bNN

fprintf('nx  nex  lam_max(S)  lam_min(S)   k(S) \n')
for nex=2:2:6;
    ney=nex;  % decomposition of Omega in nex x ney rectangles
for nx=4;     % polynomial degree in each element along y-direction
    ny=nx;
    gammax=[]; gammay=[]; % if SEM decomposition is not uniform: 
                          % they are the arrays with intefaces positions


[param]=eig_schur_2d(xa,xb,ya,yb,gam,...
          cb,nex,nx,ney,ny,gammax,gammay,param);
 
% output
fprintf('nx=%d,nex=%d,lam_max(S)=%11.4e, lam_min(S)=%11.4e, k(S)=%11.4e \n',...
    nx,nex,param(11),param(12),param(13));

end
end
