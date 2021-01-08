function plot_2dlagrange
% PLOT_2DLAGRANGE  plots 2D Lagrange polynomials.
%             Formula (1.2.55), pag. 17, CHQZ2 (tensorial product),
%             produces part of Fig. 2.13, pag. 100 CHQZ2
%
%     plot_2dlagrange   % no input, no output
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$
%

addpath ../Level_0
nx=6;
xga=-1; xgb=1.;
npdx=nx+1; 
[xgl,wx]=xwlgl(npdx); 
nx1=25;
x=linspace(-1,1,nx1)';
a=intlag_lgl(xgl, x);
ygl=zeros(npdx,1); ygl(1)=1;
y0=a*ygl;
ygl=zeros(npdx,1); ygl(3)=1;
y2=a*ygl;
[x,y]=meshgrid(linspace(-1,1,nx1));


fig=figure(...,
    'Name','2D Lagrange (vertex) polynomial',...
    'Visible','on');
graymon
z=y0*y0'; mesh(x,y,z);
xlabel('x'); ylabel('y');
view([116,32])
% print(fig,'-deps2','lagvertex2d')

fig=figure(...,
    'Name','2D Lagrange (edge) polynomial',...
    'Visible','on');
graymon
z=y0*y2'; mesh(x,y,z)
xlabel('x'); ylabel('y');
view([116,32])
%print(fig,'-deps2','lagedge2d')
% 
fig=figure(...,
    'Name','2D Lagrange (bubble) polynomial',...
    'Visible','on');
graymon
z=y2*y2'; mesh(x,y,z)
xlabel('x'); ylabel('y');
view([116,32])
%print(fig,'-deps2','lagbubble2d')



return

