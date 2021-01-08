% PLOT_MESH Script for plotting SEM mesh  on a rectangle
%
xa=-1;xb=1;   % Omega=(xa,xb) x (ya,yb)
ya=-1;yb=1;
cb='dddd';
nex=1;
    ney=nex;  % decomposition of Omega in nex x ney rectangles
nx=12;
    ny=nx;     % polynomial degree in each element along y-direction
gammax=[]; gammay=[]; 
xy=[]; 

npdx=nx+1; npdy=ny+1; ldnov=npdx*npdy; mn=ldnov; ne=nex*ney;

[x,wx]=xwlgl(npdx);
[y,wy]=xwlgl(npdy);
%
% nov construction
%
[nov]=cosnov_2d(npdx,nex,npdy,ney);
noe=nov(ldnov,ne);

% Mesh generation

[xx,yy,jacx,jacy,xy,ww,ifro]=mesh_2d(xa,xb,ya,yb,cb,nex,ney,npdx,npdy,...
nov,x,wx,y,wy,gammax,gammay);
nm=npdx*(npdy-1);
plot(xy(:,1),xy(:,2),'ko','markerfacecolor','k','markersize',4);hold on
for j=1:npdy
    plot([xy((j-1)*npdx+1,1),xy(j*npdx,1)],[y(j),y(j)],'k');
end
for i=1:npdx
    plot([x(i),x(i)],[xy(1,2),xy(nm+i,2)],'k');
end
axis equal
axis off
text(-0.14,0.16,'h','Fontsize',16)
plot([xy(85,1),xy(97,1)],[xy(85,2),xy(97,2)],'r','Linewidth',2);
hold off
print(1,'-deps2','mesh2d')
