function plot_legendre
% PLOT_LEGENDRE  plots 1D Legendre polynomials.
%             Formula  (2.3.2), pag. 75, CHQZ2,
%             produces Fig. 2.12, pag. 83 CHQZ2
%
% plot_legendre % no input, no output
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

addpath ../Level_0
l0='1'; l1='x'; 
l2=pol_legendre(1,l1,l0);
l3=pol_legendre(2,l2,l1);
l4=pol_legendre(3,l3,l2);
nx1=51;
x=linspace(-1,1,nx1)';
y=zeros(nx1,5);
y(:,1)=eval(l0); y(:,2)=eval(l1); y(:,3)=eval(l2); 
y(:,4)=eval(l3); y(:,5)=eval(l4);

for k=1:5
fig=figure(...,
    'Name','1D Legendre polynomial',...
    'Visible','on');
plot([-1.1,1.1],[0,0],'k');
hold on
plot([-1.,1.],[0,0],'k+');
plot(x,y(:,k),'k','LineWidth',3)
set(gca,'DataAspectRatioMode','manual',...
'DataAspectRatio',[1,1,1],...
'XTickLabelMode','manual','XTickLabel',[],'XTick',[],...
'YTickLabelMode','manual','YTickLabel',[],'YTick',[],...
'color','none','xcolor',[1 1 1],'ycolor',[1 1 1],...
'Box','off','FontName','Times','Fontsize',16)
if mod(k,2)==1
set(gca,'Ylimmode','manual','Ylim',[-.55,1.1])
else
set(gca,'Ylimmode','manual','Ylim',[-1.1,1.1])
end

h1=text(-1.1,-0.2,0.,'-1','FontName','Times','Fontsize',16);
h2=text(.99,-0.2,0.,'1','FontName','Times','Fontsize',16);
h3=text(.6,.8,0.,['P',num2str(k-1)],'FontName','Times','Fontsize',16);
nomefile=['figure2-legendre',num2str(k-1)];
axis off

end
return
