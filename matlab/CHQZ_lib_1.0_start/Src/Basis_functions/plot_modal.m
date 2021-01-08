function plot_modal
% PLOT_MODAL  plots 1D modal boundary-adapted polynomials 
%             Formula (2.3.31), pag. 82, CHQZ2,
%             produces part of Fig. 2.12, pag. 83 CHQZ2
%
%     plot_modal   % no input, no output
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

%
addpath ../Level_0
lm1='0';
l0='1'; l1='x'; 
l2=pol_legendre(1,l1,l0);
l3=pol_legendre(2,l2,l1);
l4=pol_legendre(3,l3,l2);
eta0='0.5*(1-x)'; 
eta1='0.5*(1+x)';
eta2=pol_modal(1,l2,l0);
eta3=pol_modal(2,l3,l1);
eta4=pol_modal(3,l4,l2);

nx1=51;
x=linspace(-1,1,nx1)';
y=zeros(nx1,5);
y(:,1)=eval(eta0); y(:,2)=eval(eta1); y(:,3)=eval(eta2); 
y(:,4)=eval(eta3); y(:,5)=eval(eta4);

for k=1:5
fig=figure(...,
    'Name','1D Modal basis polynomial',...
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
h1=text(-1.2,-0.2,0.,'-1','FontName','Times','Fontsize',16);
h2=text(.99,-0.2,0.,'1','FontName','Times','Fontsize',16);
h3=text(.7,.7,0.,['E',num2str(k-1)],'FontName','Times','Fontsize',16);
nomefile=['figure2-modal',num2str(k-1)];
axis off

end
return
