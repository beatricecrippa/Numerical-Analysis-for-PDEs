function[xx,yy,zz,jacx,jacy,jacz,xyz,ww,ifro,nov]=mesh3d(xa,xb,ya,yb,za,zb,...
nex,ney,nez,npdx,npdy,npdz,x,wx,y,wy,z,wz,gammax,gammay,gammaz);
% MESH3D   Uniform 3D Spectral element mesh on parallelepipedon
%
%           Omega=(xa,xb) x (ya,yb) x (za,zb)
%
%             V8  _________ V7
%               /|        /|
%              / |       / |
%          V5 /________ /V6|
%             |  |      |  |
%             |V4|______|__| V3
%             |  /      | /               
%             | /       |/
%             |/________/
%           V1            V2
%                      
%
%
%       __________________________
%       |      |      |     |     |
%       |  3   |  6   |  9  | 12  |      Spectral elements
%       |      |      |     |     |      ordering in a plane yz, then 
%       __________________________       planes at different x follow.
%       |      |      |     |     |
%       |  2   |  5   |  8  | 11  |
%       |      |      |     |     |
%       __________________________
%       |      |      |     |     |
%       |  1   |  4   |  7  | 10  |
%       |      |      |     |     |
%       __________________________
%
%  [xx,yy,zz,jacx,jacy,jacz,xyz,ww,ifro,nov]=mesh3d(xa,xb,ya,yb,za,zb,...
%  nex,ney,nez,npdx,npdy,npdz,nov,x,wx,y,wy,z,wz,gammax,gammay,gammaz);
%
% Input: xa= abscissa of either vertex V1 or vertex V4
%        xb= abscissa of either vertex V2 or vertex V3
%        ya= ordinate of either vertex V1 or vertex V2
%        yb= ordinate of either vertex V3 or vertex V4
%        za= ordinate of either vertex V1 
%        zb= ordinate of either vertex V5 
%        nex = number of elements along x-direction
%        ney = number of elements along y-direction
%        nez = number of elements (equally spaced) along z-direction
%        npdx = number of nodes in each element (the same in every element)
%               along x-direction
%        npdy = number of nodes in each element (the same in every element)
%               along y-direction
%        npdz = number of nodes in each element (the same in every element)
%               along z-direction
%        nov = local -global map, previously generated by cosnov_2d
%        x = npdx LGL nodes in [-1,1], previously generated by xwlgl
%        wx = npdx LGL weigths in [-1,1], previously generated by xwlgl
%        y = npdy LGL nodes in [-1,1], previously generated by xwlgl
%        wy = npdy LGL weigths in [-1,1], previously generated by xwlgl
%        z = npdz LGL nodes in [-1,1], previously generated by xwlgl
%        wz = npdz LGL weigths in [-1,1], previously generated by xwlgl
%        gammax = column or row array of length nex-1. 
%               If the deomposition of Omega is not
%               uniform, gammax is the vector of position of interfaces between
%               spectral elements along x-direction. If gammax=[], uniform
%               decomposition is used.
%        gammay = column or row array of length ney-1. 
%               If the deomposition of Omega is not
%               uniform, gammay is the vector of position of interfaces between
%               spectral elements along y-direction. If gammay=[], uniform
%               decomposition is used.
%        gammaz = column or row array of length nez-1. 
%               If the deomposition of Omega is not
%               uniform, gammay is the vector of position of interfaces between
%               spectral elements along z-direction. If gammay=[], uniform
%               decomposition is used.
%
% Output: xx = 2-indexes array of size (8,ne) 
%            xx(1:8,ie)=[x_V1_ie;x_V2_ie;x_V3_ie;x_V4_ie;...
%                        x_V5_ie;x_V6_ie;x_V7_ie;x_V8_ie]
%            (ne=nex*ney*nez is the global number of spectral elements)
%            yy(1:8,ie)=[y_V1_ie;y_V2_ie;y_V3_ie;y_V4_ie;...
%                        y_V5_ie;y_V6_ie;y_V7_ie;y_V8_ie]
%            zz(1:8,ie)=[z_V1_ie;z_V2_ie;z_V3_ie;z_V4_ie;...
%                        z_V5_ie;z_V6_ie;z_V7_ie;z_V8_ie]
%         jacx = array (length(jacx)=ne); jacx(ie)= (x_V2_ie-x_V1_ie)/2
%         jacy = array (length(jacy)=ne); jacy(ie)= (y_V3_ie-y_V1_ie)/2
%         jacz = array (length(jacz)=ne); jacz(ie)= (z_V5_ie-z_V1_ie)/2
%         xyz = column array with global mesh, length: noe=nov(npdx*npdy*npdz,ne)
%         ww = column array with global weigths, length: noe=nov(npdx*npdy*npdz,ne)
%              diag(ww) is the mass matrix associated to SEM discretization
%         ifro = column array of length noe=nov(npdx*npdy*npdz,ne): 
%            if (x_i,y_i) is internal to Omega then ifro(i)=0,
%            if (x_i,y_i) is on \partial\Omega then ifro(i)=1,
%         nov = 2-index array of local to global map,
%               size(nov)=[max(npdx*npdy*npdz),ne]
%
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

mn=npdx*npdy; ldnov=mn*npdz; ne=nex*ney*nez;
xx=zeros(8,ne); yy=zeros(8,ne); zz=zeros(8,ne);
jac=zeros(ne,1);
jacx=zeros(ne,1);
jacy=zeros(ne,1);
jacz=zeros(ne,1);
if sum(size(gammax))==0 
    Hx=(xb-xa)/nex*ones(nex);
else 
    Hx=zeros(nex,1);
    Hx(1)=gammax(1)-xa;
    Hx(2:nex-1)=gammax(2:nex-1)-gammax(1:nex-2);
    Hx(nex)=xb-gammax(nex-1);
end
if sum(size(gammay))==0 
    Hy=(yb-ya)/ney*ones(ney);
else 
    Hy=zeros(ney,1);
    Hy(1)=gammay(1)-ya;
    Hy(2:ney-1)=gammay(2:ney-1)-gammay(1:ney-2);
    Hy(ney)=yb-gammay(ney-1);
end
if sum(size(gammaz))==0 
    Hz=(zb-za)/nez*ones(nez);
else 
    Hz=zeros(nez,1);
    Hz(1)=gammaz(1)-za;
    Hz(2:nez-1)=gammaz(2:nez-1)-gammaz(1:nez-2);
    Hz(nez)=zb-gammaz(nez-1);
end

% construction of xx, yy, zz, jacx, jacy, jacz  


% first element 
      xx(1:4,1)=[xa;xa+Hx(1);xa+Hx(1);xa];
      xx(5:8,1)=xx(1:4,1);
      yy(1:2,1)=[ya,ya]; yy(3:4,1)=yy(1:2,1)+Hy(1);
      yy(5:6,1)=yy(1:2,1);yy(7:8,1)=yy(3:4,1);
      zz(1:4,1)=[za;za;za;za];
      zz(5:8,1)=[za+Hz(1);za+Hz(1);za+Hz(1);za+Hz(1)];
      jacx(1)=Hx(1)*.5;jacy(1)=Hy(1)*.5; jacz(1)=Hz(1)*.5;
      
% first row of elements along y-direction 
for iey=2:ney
      ie=iey;
      xx(1:8,ie)=xx(1:8,1);
      yy(1:2,ie)=yy(3:4,iey-1); yy(3:4,ie)=yy(1:2,ie)+Hy(iey);
      yy(5:6,ie)=yy(1:2,ie);yy(7:8,ie)=yy(3:4,ie);
      zz(1:4,ie)=zz(1:4,1); zz(5:8,ie)=zz(1:4,ie)+Hz(1);
jacx(ie)=Hx(1)*.5;
jacy(ie)=Hy(iey)*.5;
jacz(ie)=Hz(1)*.5;
end
      
% other elements on the plane xz

for iez=2:nez
for iey=1:ney
      ie=(iez-1)*ney+iey;

      xx(1:8,ie)=xx(1:8,1);
      yy(1:8,ie)=yy(1:8,iey);
      zz(1:4,ie)=zz(5:8,ie-ney); zz(5:8,ie)=zz(1:4,ie)+Hz(iez);
jacx(ie)=Hx(1)*.5;
jacy(ie)=Hy(iey)*.5;
jacz(ie)=Hz(iez)*.5;
end
end

% the whole mesh

ne2=ney*nez;
for iex=2:nex
for iez=1:nez
for iey=1:ney
      ie=(iex-1)*ne2+(iez-1)*ney+iey;
      xx([1;4;5;8],ie)=xx([2;3;6;7],ie-ne2);
      xx([2;3;6;7],ie)=xx([1;4;5;8],ie)+Hx(iex);
      yy(1:8,ie)=yy(1:8,iey);
      zz(1:8,ie)=zz(1:8,ie-ne2);
jacx(ie)=Hx(iex)*.5; 
jacy(ie)=Hy(iey)*.5; 
jacz(ie)=Hz(iez)*.5;
end
end
end

% construction of all local meshes

xyz_loc=zeros(ldnov,3,ne);
ifro_loc=zeros(ldnov,ne);
ww_loc=zeros(ldnov,ne);
for ie=1:ne
bpa=(xx(1,ie)+xx(2,ie))*.5; 
dpc=(yy(1,ie)+yy(3,ie))*.5; 
epf=(zz(1,ie)+zz(5,ie))*.5;
jac=jacx(ie)*jacy(ie)*jacz(ie);
for l=1:npdz
for j=1:npdy
for i=1:npdx
k=(l-1)*mn+(j-1)*npdx+i;
      xyz_loc(k,1,ie)=x(i)*jacx(ie)+bpa;
      xyz_loc(k,2,ie)=y(j)*jacy(ie)+dpc;
      xyz_loc(k,3,ie)=z(l)*jacz(ie)+epf;
ww_loc(k,ie)=wz(l)*wx(i)*wy(j)*jac; 
end
end

end

for l=1:npdz
for j=1:npdy
for i=1:npdx
  k=(l-1)*mn+(j-1)*npdx+i;
  if (l==1 | l==npdz | j==1 | j==npdy | i==1  | i==npdx)
  ifro_loc(k,ie)=1;
  end
end
end
end

end


% assembling global and unique mesh
xyz=zeros(ldnov*ne,3); 
nov=zeros(ldnov,ne);
xyz(1:ldnov,1:3)=xyz_loc(1:ldnov,1:3,1);
nov(1:ldnov,1)=(1:ldnov)';
ifro(1:ldnov)=ifro_loc(1:ldnov,1);
epsi=1.e-14;
last=ldnov;
for ie=2:ne
for i=1:ldnov
punto_new(1:3)=xyz_loc(i,1:3,ie);

tt=0;j=1;
while (j<=last & tt==0)
if (ifro(j)~=0)

dd=abs(punto_new-xyz(j,1:3));
if (dd<=epsi)
nov(i,ie)=j;
tt=1;
end
end
j=j+1;
end
if (nov(i,ie)==0)
last=last+1;
xyz(last,1:3)=punto_new(1:3);
ifro(last)=ifro_loc(i,ie);
nov(i,ie)=last;
end

end
end

ifro=ifro';
noe=last;
xyz=xyz(1:noe,1:3);
ww=zeros(noe,1);

for ie=1:ne
    ww(nov(1:ldnov,ie))=ww(nov(1:ldnov,ie))+ww_loc(1:ldnov,ie);
end
clear ww_loc xyz_loc ifro_loc

for i=1:noe
      if (ifro(i)~=0)
      if((abs(xyz(i,1)-xa) > eps & abs(xyz(i,1)-xb) > eps)...
     & (abs(xyz(i,2)-ya) > eps & abs(xyz(i,2)-yb) > eps)...
     & (abs(xyz(i,3)-za) > eps & abs(xyz(i,3)-zb) > eps))
      ifro(i)=0;
      end
      end
end


