% FIG4_14  Script to produce Fig 4.14,  pag. 206 CHQZ2
%
% Legendre collocation first-derivative eigenvalues and spectral condition
% number.
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$


N=(8:8:64);
    nu=1;
% left
dmin0=[];dmax0=[]; k20=[];
dmin1=[];dmax1=[]; k21=[];
dmin2=[];dmax2=[]; k22=[];
for nx=N
    pbl=0;
    [d,A]=lgl_eig(nx,nu,pbl);
    d=abs(d);
    dmin0=[dmin0;min(d)];
    dmax0=[dmax0;max(d)];
    k20=[k20;cond(A)];
    pbl=1;
    [d,A]=lgl_eig(nx,nu,pbl);
    d=abs(d);
    dmin1=[dmin1;min(d)];
    dmax1=[dmax1;max(d)];
    k21=[k21;cond(A)];
    pbl=2;
    [d,A]=lgl_eig(nx,nu,pbl);
    d=abs(d);
    dmin2=[dmin2;min(d)];
    dmax2=[dmax2;max(d)];
    k22=[k22;cond(A)];
end
fig=figure(...,
    'Name','Fig. 4.14 left',...
    'Visible','on')
loglog(N,dmin0,'k-o',N,dmax0,'k-+')
hold on
loglog(N,dmin1,'k->',N,dmax1,'k-s')
loglog(N,dmin2,'k-*',N,dmax2,'k-d')
set(gca,'Xlimmode','manual','Xlim',[4,130],...
    'Ylimmode','manual','Ylim',[1.e-2,1.e6],...
    'Xgrid','on','XminorGrid','off','Ygrid','on','YminorGrid','off',...
'LineWidth',1,...
'FontName','Times','Fontsize',16,'FontWeight','normal')

fig=figure(...,
    'Name','Fig. 4.14 right',...
    'Visible','on')
loglog(N,k20,'k-o')
hold on
loglog(N,k21,'k->')
loglog(N,k22,'k-*')
set(gca,'Xlimmode','manual','Xlim',[4,130],...
    'Ylimmode','manual','Ylim',[10,1.e5],...
    'Xgrid','on','XminorGrid','off','Ygrid','on','YminorGrid','off',...
'LineWidth',1,...
'FontName','Times','Fontsize',16,'FontWeight','normal')
