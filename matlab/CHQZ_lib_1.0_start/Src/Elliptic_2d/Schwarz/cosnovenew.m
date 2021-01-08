function [nove,nvle]=cosnovenew(nx,nex,ny,ney,nov,ifro,nlevel);
%  COSNOVENEW   Construction of restriction maps for extended elements 
%
%      [nove,nvle]=cosnovenew(npdx,nex,npdy,ney,nov,ifro,nlevel);
%
% Input: nx = polynomial degree in each spectral element along x-direction
%        nex = number of spectral elements along x-direction
%        ny = polynomial degree in each spectral element along y-direction
%        ney = number of spectral elements along y-direction
%        nov = 2-index array of local to global map, 
%               size(nov)=[max(npdx*npdy),ne]
%        ifro = column array of length noe=nov(npdx*npdy,ne): 
%            if (x_i,y_i) is internal to Omega then ifro(i)=0,
%            if (x_i,y_i) is on \partial\Omega then ifro(i)=1,
%        nlevel = number of added layers for extending spectral elements
%                   inside additive Schwarz preconditioner.
%                   If param(2)==1, the preconditioner is P^{m}_{as,H} 
%                                    (minimum overlap), pag. 377 CHQZ3)
%                   If param(2)==2, the preconditioner is P^{s}_{as,H} 
%                                    (small overlap), pag. 377 CHQZ3)
%                   param(2) is a positive integer less than min(nx,ny)
%
% Output: nove = 2-index array of "extended local" to global map, 
%         nvle = 2-index array:
%                 nvle (:,1)=number of nodes of the extendend elements
%                 nvle (:,2)=number of Q1 elements of the extendend elements
%                            along x-direction
%                 nvle (:,3)=number of Q1 elements of the extendend elements
%                            along y-direction
%
%
% References: CHQZ3 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Evolution to Complex Geometries 
%                     and Applications to Fluid DynamicsSpectral Methods"
%                    Springer Verlag, Berlin Heidelberg New York, 2007.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

[ldnov,ne]=size(nov);

npdx=nx+1; npdy=ny+1;mn=npdx*npdy;
nm1=mn-nx;
mne=(npdx+nlevel*2)*(npdy+nlevel*2);
nove=zeros(mne,ne);
nvle=zeros(ne,3);
for ie=1:ne
nove(1:mn,ie)=nov(1:mn,ie);
nvle(ie,2)=nx; nvle(ie,3)=ny;
ultimo=mn;
vertex=zeros(4,1);

% vertex 1

ig=nov(1,ie);
if(ifro(ig)<=0)
vertex(1)=1;
% The vertex is inside Omega, I expand nov
% by adding the last  "nlevel" rows of the domain down
for ib=1:nlevel
    nove(ultimo+1:ultimo+npdx,ie)=nov(npdx*(ny-ib)+1:npdx*(ny-ib+1),ie-1);
ultimo=ultimo+npdx; nvle(ie,3)=nvle(ie,3)+1;
end

% by adding the last  "nlevel" columns of the left domain
for ib=1:nlevel
nove(ultimo+1:ultimo+npdy,ie)=nov(nx-ib+1:npdx:mn-ib,ie-ney);
ultimo=ultimo+npdy; nvle(ie,2)=nvle(ie,2)+1;
end
% by adding nodes of the element symmetric to "ie" with respect to the vertex
for ib=1:nlevel
    nove(ultimo+1:ultimo+nlevel,ie)=nov(mn-npdx*ib-nlevel:mn-npdx*ib-1,ie-ney-1);
ultimo=ultimo+nlevel;
end
end

% vertex 2

ig=nov(npdx,ie);
if(ifro(ig)<=0)
vertex(2)=1;
% The vertex is inside Omega, I expand nov
if vertex(1)==0
% by adding the last  "nlevel" rows of the  domain down
for ib=1:nlevel
    nove(ultimo+1:ultimo+npdx,ie)=nov(npdx*(ny-ib)+1:npdx*(ny-ib+1),ie-1);
ultimo=ultimo+npdx;nvle(ie,3)=nvle(ie,3)+1;
end
end
% by adding the first  "nlevel" columns of the right domain
for ib=1:nlevel
nove(ultimo+1:ultimo+npdy,ie)=nov(ib+1:npdx:nm1+ib,ie+ney);
ultimo=ultimo+npdy; nvle(ie,2)=nvle(ie,2)+1;
end
% by adding nodes of the element symmetric to "ie" with respect to the vertex
 for ib=1:nlevel
     nove(ultimo+1:ultimo+nlevel,ie)=nov(mn-npdx*(ib+1)+2:mn-npdx*(ib+1)+1+nlevel,ie+ney-1);
 ultimo=ultimo+nlevel;
end
end

% vertex 3

ig=nov(mn,ie);
if(ifro(ig)<=0)
vertex(3)=1;
% The vertex is inside Omega, I expand nov
if vertex(2)==0
% by adding the first  "nlevel" columns of the right domain
for ib=1:nlevel
    nove(ultimo+1:ultimo+npdy,ie)=nov(ib+1:npdx:nm1+ib,ie+ney);
ultimo=ultimo+npdy; nvle(ie,2)=nvle(ie,2)+1;
end
end
% by adding the last  "nlevel" rows of the  domain up
for ib=1:nlevel
    nove(ultimo+1:ultimo+npdx,ie)=nov(npdx*ib+1:npdx*(ib+1),ie+1);
ultimo=ultimo+npdx; nvle(ie,3)=nvle(ie,3)+1;
end
% by adding nodes of the element symmetric to "ie" with respect to the vertex
for ib=1:nlevel
    nove(ultimo+1:ultimo+nlevel,ie)=nov(npdx*ib+2:npdx*ib+1+nlevel,ie+ney+1);
ultimo=ultimo+nlevel;
end
end

% vertex 4

ig=nov(nm1,ie);
if(ifro(ig)<=0)
vertex(4)=1;
% The vertex is inside Omega, I expand nov
if vertex(1)==0
% by adding the last  "nlevel" columns of the left domain
for ib=1:nlevel
    nove(ultimo+1:ultimo+npdy,ie)=nov(nx-ib+1:npdx:mn-ib,ie-ney);
ultimo=ultimo+npdy; nvle(ie,2)=nvle(ie,2)+1;
end
end
if vertex(3)==0
% by adding the last  "nlevel" rows of the  domain up
for ib=1:nlevel
    nove(ultimo+1:ultimo+npdx,ie)=nov(npdx*ib+1:npdx*(ib+1),ie+1);
ultimo=ultimo+npdx; nvle(ie,3)=nvle(ie,3)+1;
end
end
% by adding nodes of the element symmetric to "ie" with respect to the vertex
for ib=1:nlevel
    nove(ultimo+1:ultimo+nlevel,ie)=nov(npdx*(ib+1)-nlevel:npdx*(ib+1)-1,ie-ney+1);
ultimo=ultimo+nlevel;
end
end
nvle(ie,1)=ultimo;

end  % end loop ie
