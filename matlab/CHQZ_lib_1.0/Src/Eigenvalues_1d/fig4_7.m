% FIG4_7  Script to produce Fig 4.7, pag. 199 CHQZ2. 
%
% Extreme eigenvalues of 
% Legendre G-NI stiffness matrices for the 2nd order derivative operator
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$


dmin=[];dmax=[];
N=(8:8:96);
    nu=1;
    pbl=31;
% left
for nx=N
    [d,A]=lgl_eig(nx,nu,pbl);
    dmin=[dmin;min(d)];
    dmax=[dmax;max(d)];
end
fig=figure(...,
    'Name','Fig. 4.7 left',...
    'Visible','on')
loglog(N,dmin,'k-.',N,dmax,'k-')
set(gca,'Xlimmode','manual','Xlim',[5,100],...
    'Ylimmode','manual','Ylim',[1.e-4,1.e6],...
    'Xgrid','on','XminorGrid','off','Ygrid','on','YminorGrid','off',...
'LineWidth',1,...
'FontName','Times','Fontsize',16,'FontWeight','normal')

% right
dmin1=[];dmax1=[];
pbl=41;
for nx=N
    [d,A]=lgl_eig(nx,nu,pbl);
temp=[];
for i=1:length(A);
if abs(d(i))>1.d-12
temp=[temp;d(i)];
end
end
dmin1=[dmin1;min(abs(temp))];
dmax1=[dmax1;max(d)];
end
fig=figure(...,
    'Name','Fig. 4.7 right',...
    'Visible','on')
loglog(N,dmin1,'k-.',N,dmax1,'k-')
set(gca,'Xlimmode','manual','Xlim',[5,100],...
    'Ylimmode','manual','Ylim',[1.e-4,1.e6],...
    'Xgrid','on','XminorGrid','off','Ygrid','on','YminorGrid','off',...
'LineWidth',1,...
'FontName','Times','Fontsize',16,'FontWeight','normal')

