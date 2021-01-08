% CALL_ADR_1D Script for producing  data of figure 1.4, pag. 21, CHQZ2
%
% to produce data of figure 1.4, pag. 21, CHQZ2:  
%>> p=1; call_adr_1d;
%>> p=2; call_adr_1d;
%>> p=3; call_adr_1d;
%>> p=4; call_adr_1d;

xa=-1;xb=1; 

% sets functions 

[uex,uexx,ff,nu,beta,gam]=setfun_adr_1d;

if p==1
nne=[10,20,40,80,160,320];
nnx=1;
elseif p==2
nne=[10,20,40,80,160];
nnx=2;
elseif p==3
nne=[10,20,40,80,108];
nnx=3;
elseif p==4
nne=1; nnx=4:8:40;
end

param=zeros(10,1);
nomefile=['dt2',num2str(p),'.dat'];
cb='dn'; % c.b.: the first character specifies boundary condition in xa
        % the second character specifies boundary condition in xb
fid=fopen(nomefile,'w');
fprintf(fid,'1D Advection-diffusion-reaction problem, FEM/SM\n');
fprintf(fid,'p  noe   e_inf  e_l2   e_h1    der   \n');

for ne=nne; %[10,20,40,80,160] % number of spectral elements
for nx=nnx % spectral polynomial degree (the same in every element)
    H=(xb-xa)/ne;
param(1:6)=[1,0,nx*2,1,0,64]; % see help adr_1d
[xy,un,err_inf,err_l2,err_h1,der]=adr_1d(xa,xb,nu,beta,gam,uex,uexx,ff,cb,ne,p,nx,param);
noe=length(xy);
fprintf('p=%d, noe=%d, e_inf=%8.2e, e_l2=%8.2e, e_h1=%8.2e, der=%8.2e \n',...
    p,noe,err_inf,err_l2,err_h1,der);
fprintf(fid,'%d   %d   %8.2e   %8.2e    %8.2e  %8.2e\n',...
    p, noe,err_inf,err_l2,err_h1,der);
end
end
fclose(fid);

