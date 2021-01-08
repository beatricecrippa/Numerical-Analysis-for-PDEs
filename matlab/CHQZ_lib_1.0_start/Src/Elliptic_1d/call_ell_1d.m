% CALL_ELL_1D Script to produce results of fig. 4.17, pag. 207, CHQZ2
%
% -nu d^2 u/dx^2 +du/dx=1  in (-1,1)
%
%   u(-1)=u(1)=0
%
%
%  Generalized Numerical integration, one element
%

xa=-1;xb=1;  
ub=[0,0];
nu=1.d-2; beta=1;  gam=0; ff=@(x)[1];
ne=1; cb='dd';
nnx=[12,20,72]; % for nu=1.d-2

fig=figure(...,
'Name', 'Advection-diffusion solution nu=1.d-2',...
'Visible','on')
color=['k-.';'k--';'k- '];
hold on
knx=0;
for nx=nnx
[xy,un]=ell_1d(xa,xb,nu,beta,gam,ff,cb,ub,ne,nx);
knx=knx+1;
% evaluation of numerical solution at a finer grid to plot it
nx1=100;
x_int=xwlgl(nx1);
[u_int]=legendre_tr_eval(xy,un,x_int);
% evaluation of the polynomial at x_int (finer mesh with respect to x)
if(knx<=3)
h(knx)=plot(x_int,u_int,color(knx,:));
else
h(knx)=plot(xx,u1,color(knx,:));
end
end
hl1=legend('N=12','N=20','N=72');
set(hl1,'Position',[0.145,0.621,0.25,0.209]);
set(gca,'XlimMode','manual','Xlim',[-1,1.001],...
'YlimMode','manual','Ylim',[0,3])
set(hl1,'FontName','Times','Fontsize',16);
set(h(1),'Linewidth',1);
set(h(2),'Linewidth',1);
set(h(3),'Linewidth',2);
set(gca,'PlotBoxAspectRatio',[3 2 1],...
'Xgrid','on','XminorGrid','off','Ygrid','on',...
'YminorGrid','off','LineWidth',1,...
'FontName','Times','Fontsize',16)
xlabel('x','FontName','Times','Fontsize',16)
hold off


nu=1.d-3; beta=1;  gam=0; ff=@(x)[1];
ne=1;
nnx=[48,72,104]; % for nu=1.d-2

fig=figure(...,
'Name', 'Advection-diffusion solution nu=1.d-3',...
'Visible','on')
color=['k-.';'k--';'k- '];
hold on
knx=0;
for nx=nnx
[xy,un]=ell_1d(xa,xb,nu,beta,gam,ff,cb,ub,ne,nx);
knx=knx+1;
% evaluation of numerical solution at a finer grid to plot it
nx1=300;
x_int=xwlgl(nx1);
[u_int]=legendre_tr_eval(xy,un,x_int);
% evaluation of the polynomial at x_int (finer mesh with respect to x)
if(knx<=3)
h(knx)=plot(x_int,u_int,color(knx,:));
else
h(knx)=plot(xx,u1,color(knx,:));
end
end
hl1=legend('N=48','N=72','N=104');
set(hl1,'Position',[0.145,0.621,0.25,0.209]);
set(gca,'XlimMode','manual','Xlim',[-1,1.001],...
'YlimMode','manual','Ylim',[0,3])
set(hl1,'FontName','Times','Fontsize',16);
set(h(1),'Linewidth',1);
set(h(2),'Linewidth',1);
set(h(3),'Linewidth',2);
set(gca,'PlotBoxAspectRatio',[3 2 1],...
'Xgrid','on','XminorGrid','off','Ygrid','on',...
'YminorGrid','off','LineWidth',1,...
'FontName','Times','Fontsize',16)
xlabel('x','FontName','Times','Fontsize',16)
hold off

return
