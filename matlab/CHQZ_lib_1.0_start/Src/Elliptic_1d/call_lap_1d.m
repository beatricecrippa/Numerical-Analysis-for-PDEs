% CALL_LAP_1D Script for calling lap_1d

xa=-1;xb=1; 

% set functions 

[uex,uexx,ff,nu,gam]=setfun_lap_1d;

nomefile=['out.dat'];
cb='nd'; % c.b.: the first character specifies boundary condition in xa
        % the second character specifies boundary condition in xb
fid=fopen(nomefile,'w');
fprintf(fid,'1D Laplace problem, SEM\n');
fprintf(fid,'nx, ne   e_inf  e_l2   e_h1 \n');
param=zeros(10,1);

nnx=4:4:24;
for ne=10; %[10,20,40,80,160] % number of spectral elements
for nx=nnx % spectral polynomial degree (the same in every element)
    H=(xb-xa)/ne;
param(1:6)=[1,1,nx*2,1,1,64]; % see help adr_1d
[xy,un,err_inf,err_l2,err_h1]=lap_1d(xa,xb,nu,gam,uex,uexx,ff,cb,ne,nx,param);
noe=length(xy);
fprintf('nx=%d, ne=%d, e_inf=%8.2e, e_l2=%8.2e, e_h1=%8.2e \n',...
    nx,ne,err_inf,err_l2,err_h1);
fprintf(fid,'%d   %d   %8.2e   %8.2e    %8.2e \n',...
    nx, ne,err_inf,err_l2,err_h1);
end
end
fclose(fid);

