% CALL_ELLPRECOFEM_1D Script to set data for calling ellprecofem_1d. 
%
% See help of ellprecofem_1d.m for a description of inputs

xa=-1;xb=1;
param=zeros(20,1);
param(1)=1;  % preconditioner
param(2)=11; % solver or eigenvalues computation
param(3)=1.d-14; % tolerance
param(4)=400;    % maximum number of iterations
[uex,uexx,ff,nu,beta,gam]=setfun_ell_1d;
ne=1;
cb='dd';  ub=[uex(xa),uex(xb)];
% cb='nn';  ub=[uexx(xa),uexx(xb)];
% cb='dn';  ub=[uex(xa),uexx(xb)];
% cb='nd';  ub=[uexx(xa),uex(xb)];
% ub=[0,0];
% fprintf('nx       ne   err_h1            iter           res\n');
if(param(2)>10)
    fprintf('nx   ne   kapit\n');
else
    fprintf('nx   ne   err_h1      iter       res\n');
end
for nx=16:16:48;
[xy,un,A,M,AFE,MFE,MFEd,d,kappa,param]=ellprecofem_1d(xa,xb,nu,beta,gam,ff,cb,...
    ub,ne,nx,param);
if(param(2)>10)
    fprintf('%d   %d   %13.6e \n',nx,ne,kappa)
elseif(param(2)<=10)
    % errors computation
param1=zeros(10,1);
param1(1:6)=[1,0,nx*2,1,1,64]; % see help ellprecofem_1d
[err_inf,err_h1,err_l2]=errors_1d(nx,ne,xa,xb,un,uex,uexx,param1);
fprintf('%d   %d   %13.6e   %d     %13.6e \n',nx,ne,err_h1,param(11),param(12))
end    
end

return
