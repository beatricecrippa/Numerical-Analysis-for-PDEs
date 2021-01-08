% CALL_STAT_HYP. Script to set input data and call  stat_scalar_hyp
%
% to set input data: help stat_scalar_hyp
% these data are useful to reproduce figure 3.7, pag. 150 CHQZ2.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

xa=-1;xb=1;
beta=1.5;
param(1)=1; param(2)=0.75;

%
%  set functions
%
uex=@(x)[sin(6*x)];
ul=sin(-6);
f=@(x)[9*cos(6*x)];
fid=fopen('error_stat','w');

% strong, penalty and GNI
for method=1:3
param(1)=method;
for nx=8:2:20;


[x,u,err]=stat_scalar_hyp(xa,xb,beta,f,uex,ul,nx,param);
fprintf('method=%d, nx=%d, err=%13.6e\n',method,nx,err);
fprintf(fid,'method=%d, nx=%d, err=%13.6e\n',method,nx,err);
end
end

fclose(fid);
