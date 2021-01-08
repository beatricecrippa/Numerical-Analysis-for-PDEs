% CALL_EIG_SCHWARZ_2D  Script file for pre and post processing eig_schwarz_2d.m
%
%

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

xa=-1;xb=1;   % Omega=(xa,xb) x (ya,yb)
ya=-1;yb=1;
cb='dddd'; % eig_schwarz_2d works only if cb='dddd';
    gammax=[]; gammay=[]; % if SEM decomposition is not uniform: 
                          % they are the arrays with intefaces positions
    param=zeros(20,1);  
    param(1)=2;    % 1:P=I, 2: P=P_as
    param(2)=2;    % number of levels for overlapping elements
fprintf('nx  nex  lam_max(S)  lam_min(S)   k(S) \n')
for nx=8  % polynomial degree in each element along x-direction
    ny=nx;     % polynomial degree in each element along y-direction
for nex=3:10;
    ney=nex;  % decomposition of Omega in nex x ney rectangles

   % call eig_schwarz_2d

[param]=eig_schwarz_2d(xa,xb,ya,yb,gam,...
          cb,nex,nx,ney,ny,gammax,gammay,param);
 
% output
fprintf('nx=%d,nex=%d,l_max(P^(-1)A)=%11.4e, l_min(P^(-1)A)=%11.4e, k(P^(-1)A)=%11.4e \n',...
    nx,nex,param(11),param(12),param(13));

end
end
