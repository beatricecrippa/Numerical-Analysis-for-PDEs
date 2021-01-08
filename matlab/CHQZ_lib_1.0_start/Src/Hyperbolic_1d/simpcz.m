function [int]=simpcx(xx,yy)
% SIMPCX. Composite Simpson Quadrature Formula 
% function [int]=simpc(xx,yy)
%
%
% Input:
% xx    equispaced abscissas, size(xx) has to be odd
% yy    ordinates 
%
% Output:
% int  computed integral 
%
% Autori:  pg
h=xx(3)-xx(1);
n=length(yy);
int=0;
m=fix(n/2);
for i=1:2:n-2
int=int+yy(i)+yy(i+2)+4*yy(i+1);
end
int=int*h/6;
