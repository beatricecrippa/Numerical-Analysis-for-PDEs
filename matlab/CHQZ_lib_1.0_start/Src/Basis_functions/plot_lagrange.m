function plot_lagrange
% PLOT_LAGRANGE  plots 1D Lagrange polynomials.
%             Formula (1.2.55), pag. 17, CHQZ2,
%             produces part of Fig. 2.12, pag. 83 CHQZ2
%
%     plot_lagrange   % no input, no output
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$
%

%
addpath ../Level_0
nx=4;
xga=-1; xgb=1.;
npdx=nx+1; 
[xgl,wx]=xwlgl(npdx); 
nx1=51;
y=zeros(nx1,5);
x=linspace(-1,1,nx1)';
a=intlag_lgl(xgl, x);

for i=1:npdx
ygl=zeros(npdx,1); ygl(i)=1;
y(:,i)=a*ygl;

fig=figure(...,
    'Name','1D Lagrange polynomial',...
    'Visible','on');
plot([-1.2,1.2],[0,0],'k'); hold on
plot(xgl,zeros(npdx,1),'k+');
plot(x,y(:,i),'k','LineWidth',3)
set(gca,'DataAspectRatioMode','manual',...
'DataAspectRatio',[1,1,1],...
'XTickLabelMode','manual','XTickLabel',[],'XTick',[],...
'YTickLabelMode','manual','YTickLabel',[],'YTick',[],...
'color','none','xcolor',[1 1 1],'ycolor',[1 1 1],...
'Box','off','FontName','Times','Fontsize',16)
if mod(i,2)==1
set(gca,'Ylimmode','manual','Ylim',[-.55,1.1])
else
set(gca,'Ylimmode','manual','Ylim',[-1.1,1.1])
end
h1=text(-1.2,-0.2,0.,'-1','FontName','Times','Fontsize',16);
h2=text(.99,-0.2,0.,'1','FontName','Times','Fontsize',16);
h3=text(.57,.8,0.,['L',num2str(i-1)],'FontName','Times','Fontsize',16);
axis off
nomefile=['figure2-lagrange',num2str(i-1)];
%print(fig,'-depsc2',nomefile)

end
return
