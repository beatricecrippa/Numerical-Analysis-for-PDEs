function [A,f]=patch_se(A,f,ifro,nov,dx,jacx,dy,jacy);
% PATCH_SE Imposes strong continuity of normal derivatives across interfaces
%
% called by lap_2d
%

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

[ldnov,ne]=size(nov);
noe=length(f); 
npdx=length(dx);
npdy=length(dy);


%  clear the rows of A related to interface unknowns

for i=1:noe
if (ifro(i)==-1)
A(i,:)=zeros(1,noe);
f(i)=0;
end
end

% derivative
mn=npdx*npdy;
for ie=1:ne
Al=patch_sp(dx,jacx(ie),dy,jacy(ie));
for j=1:npdy
for i=1:npdx
k=(j-1)*npdx+i;
ii=nov(k,ie);
if(ifro(ii)==-1)
A(ii,nov(1:mn,ie))=A(ii,nov(1:mn,ie))+Al(k,:);
end
end
end
end
return
