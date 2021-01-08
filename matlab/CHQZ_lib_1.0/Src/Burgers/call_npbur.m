% CALL_NPBUR Script to set input and to call nonperiodic_burgers.m
%            Read help of nonperiodic_burgers.m to set input data
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

Err=[];
xa=-1;xb=1;
visc=0.01;
tscheme=4;bc=1;
t0=0;T=1;
deltat=1.e-2;
NX=32;
%NX=[16,32,64,128];
for nx=NX
[u,err]=nonperiodic_burgers(xa,xb,t0,T,visc,nx,deltat,tscheme,bc);
Err=[Err;err];
end
Err
