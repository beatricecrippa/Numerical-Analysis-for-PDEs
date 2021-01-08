% CALL_HYP. Script to set input data and call  scalar_hyp and stag_scalar_hyp
%
% to set input data: help scalar_hyp, help stag_scalar_hyp
% these data are useful to reproduce figure 3.5, pag. 150 CHQZ2.
%

%   Written by Paola Gervasio
%   $Date: 2007/04/01$


xa=-1;xb=2;
beta=1.5;
param(1)=1; param(2)=0.75;
t0=0;T=1.d-3;
deltat=1.d-4;
%
%  set functions
%
uex=@(x,t)[sin(2*x-3*t)];
u0=@(x)[sin(2*x)]; ul=@(t)[sin(-2-3*t)];

fid=fopen('errorit4','w');

% strong, penalty and GNI
for method=1:0
param(1)=method;
for nx=8:2:20;


[x,u,err,Psi,Phi]=scalar_hyp(xa,xb,t0,T,beta,uex,u0,ul,nx,deltat,param);
fprintf('method=%d, nx=%d, err=%13.6e\n',method,nx,err);
fprintf(fid,'method=%d, nx=%d, err=%13.6e\n',method,nx,err);
end
end

% staggered grids
method=4;
for nx=8:2:20;

[x,u,err,Psi,Phi]=stag_scalar_hyp(xa,xb,t0,T,beta,uex,u0,ul,nx,deltat);
fprintf('method=%d, nx=%d, err=%13.6e\n',method,nx,err);
fprintf(fid,'method=%d, nx=%d, err=%13.6e\n',method,nx,err);
end
fclose(fid);
