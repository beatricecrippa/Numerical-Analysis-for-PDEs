function A=patch_sp(dx,jacx,dy,jacy);
% PATCH_SP  Computes normal derivatives on the boundary of the spectral element
%

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

npdx=length(dx);
npdy=length(dy);

nm=npdx*(npdy-1);
mn=npdx*npdy;
A=sparse(mn,mn);
% side 1
for i=1:npdx
ki=i;
for j=1:npdy
kj=(j-1)*npdx+i;
A(ki,kj)=-dy(1,j)/jacy;
end
end

% side 2
for j=1:npdy
kj=j*npdx;
for i=1:npdx
ki=(j-1)*npdx+i;
A(kj,ki)=A(kj,ki)+dx(npdx,i)/jacx;
end
end

% side 3
for i=1:npdx
ki=nm+i;
for j=1:npdy
kj=(j-1)*npdx+i;
A(ki,kj)=A(ki,kj)+dy(npdy,j)/jacy;
end
end

% side 4
for j=1:npdy
kj=(j-1)*npdx+1;
for i=1:npdx
ki=(j-1)*npdx+i;
A(kj,ki)=A(kj,ki)-dx(1,i)/jacx;
end
end

