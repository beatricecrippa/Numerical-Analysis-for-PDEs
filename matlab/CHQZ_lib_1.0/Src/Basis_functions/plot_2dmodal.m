function plot_2dmodal
% PLOT_2DMODAL  plots 2D modal boundary adapted polynomials.
%             Formula (2.3.31), pag. 82, CHQZ2 (tensorial product),
%             produces part of Fig. 2.13, pag. 100 CHQZ2
%
%     plot_2dmodal   % no input, no output
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$
%
addpath ../Level_0
nx=6;
lm1='0';
l0='1'; l1='x';
l2=pol_legendre(1,l1,l0);
l3=pol_legendre(2,l2,l1);
l4=pol_legendre(3,l3,l2);
l5=pol_legendre(3,l4,l3);
l6=pol_legendre(3,l5,l4);
eta0='0.5*(1-x)';
eta1=pol_modal(1,l2,l0);
eta2=pol_modal(2,l3,l1);
eta3=pol_modal(3,l4,l2);
eta4=pol_modal(4,l5,l3);
eta5=pol_modal(5,l6,l4);
eta6='0.5*(1+x)';

nx1=25;
x=linspace(-1,1,nx1)';
y0=eval(eta0); y2=eval(eta2);

[x,y]=meshgrid(linspace(-1,1,nx1));
fig=figure(...,
    'Name','2D Modal (vertex) polynomial',...
    'Visible','on');

graymon
z=y0*y0'; mesh(x,y,z)
xlabel('x'); ylabel('y')
view([116,32])
print(fig,'-deps2','modvertex2d')


fig=figure(...,
    'Name','2D Modal (edge) polynomial',...
    'Visible','on');
graymon
z=y0*y2'; mesh(x,y,z)
xlabel('x'); ylabel('y')
view([116,32])
print(2,'-deps2','modedge2d')


fig=figure(...,
    'Name','2D Modal (bubble) polynomial',...
    'Visible','on');
graymon
z=y2*y2'; mesh(x,y,z)
xlabel('x'); ylabel('y')
view([116,32])
print(3,'-deps2','modbubble2d')



return

