% CALL_EIG_SCHUR_2D_FILE Script for pre and post processing eig_schur_2d.m
%

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

xa=-1;xb=1;   % Omega=(xa,xb) x (ya,yb)
ya=-1;yb=1;
gam=0;
cb='dddd';     %  schur_2d works only if cb='dddd';

%dependence on N

fid=fopen('conds_schurPN_new','w');
fprintf(fid,'cond(P^-1 Sigma) per -Lap. nex=ney=2; Omega=(-1,1)^2\n');
fprintf(fid,'N,   P=I   P=NN    P=bNN   \n');
% 
fprintf('nx  nex  lam_max(S)  lam_min(S)   k(S) \n')
for nex=2; 
    ney=nex;  % decomposition of Omega in nex x ney rectangles
for nx=4:4:64;     % polynomial degree in each element along y-direction
    ny=nx;
    gammax=[]; gammay=[]; % if SEM decomposition is not uniform: 
                         % they are the arrays with intefaces positions
cc=[nx];
param=zeros(20,1); param(1)=1;    % 1:P=I, 2: P=NN, 3: P=bNN
[param]=eig_schur_2d(xa,xb,ya,yb,gam,...
          cb,nex,nx,ney,ny,gammax,gammay,param);
cc=[cc,param(13)];
param=zeros(20,1); param(1)=2;    % 1:P=I, 2: P=NN, 3: P=bNN
[param]=eig_schur_2d(xa,xb,ya,yb,gam,...
          cb,nex,nx,ney,ny,gammax,gammay,param);
cc=[cc,param(13)];
param=zeros(20,1); param(1)=3;    % 1:P=I, 2: P=NN, 3: P=bNN
[param]=eig_schur_2d(xa,xb,ya,yb,gam,...
          cb,nex,nx,ney,ny,gammax,gammay,param);
cc=[cc,param(13)];
% output
% fprintf('nx=%d,nex=%d,lam_max(S)=%11.4e, lam_min(S)=%11.4e, k(S)=%11.4e \n',...
%     nx,nex,param(11),param(12),param(13));
fprintf(fid,'%d   %9.3e  %9.3e   %9.3e  \n',cc);
fprintf('%d   %9.3e  %9.3e   %9.3e  \n',cc);
end
end
fclose(fid);

%dependence on H

fid=fopen('conds_schurPH_new','w');
fprintf(fid,'cond(P^-1 Sigma) per -Lap. nx=ny=4; Omega=(-1,1)^2\n');
fprintf(fid,'H,   P=I   P=NN    P=bNN   \n');

fprintf('nx  nex  lam_max(S)  lam_min(S)   k(S) \n')
for nex=2:2:22;
    ney=nex;  % decomposition of Omega in nex x ney rectangles
for nx=4;     % polynomial degree in each element along y-direction
    ny=nx;
    gammax=[]; gammay=[]; % if SEM decomposition is not uniform: 
                         % they are the arrays with intefaces positions
cc=[2/nex];
param=zeros(20,1); param(1)=1;    % 1:P=I, 2: P=NN, 3: P=bNN
[param]=eig_schur_2d(xa,xb,ya,yb,gam,...
          cb,nex,nx,ney,ny,gammax,gammay,param);
cc=[cc,param(13)];
param=zeros(20,1); param(1)=2;    % 1:P=I, 2: P=NN, 3: P=bNN
[param]=eig_schur_2d(xa,xb,ya,yb,gam,...
          cb,nex,nx,ney,ny,gammax,gammay,param);
cc=[cc,param(13)];
param=zeros(20,1); param(1)=3;    % 1:P=I, 2: P=NN, 3: P=bNN
[param]=eig_schur_2d(xa,xb,ya,yb,gam,...
          cb,nex,nx,ney,ny,gammax,gammay,param);
cc=[cc,param(13)];
% output
% fprintf('nx=%d,nex=%d,lam_max(S)=%11.4e, lam_min(S)=%11.4e, k(S)=%11.4e \n',...
%     nx,nex,param(11),param(12),param(13));
fprintf(fid,'%d   %9.3e  %9.3e   %9.3e  \n',cc);
fprintf('%d   %9.3e  %9.3e   %9.3e  \n',cc);
end
end
fclose(fid);
