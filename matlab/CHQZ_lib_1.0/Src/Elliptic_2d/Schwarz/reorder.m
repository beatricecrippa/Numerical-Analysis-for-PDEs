function [p]=reorder(xyl,nex,ney);
% REORDER Reordering of the array of nodes of extended element.
% Lexicographic ordering is performed
%
% [p]=reorder(xyl,nex,ney);
%
% Input : xyl = nodes of the extended element
%         nex = number of Q1 elements of the extended macro element, 
%               alobng x-direction
%         ney = number of Q1 elements of the extended macro element, 
%               alobng y-direction
%  
% Output : p = reordering array of nodes

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

npx=nex+1;npy=ney+1;
noe=npx*npy;

ctx=sort(xyl(:,1)); 
ctx=ctx(1:npy:noe);
cty=sort(xyl(:,2)); 
cty=cty(1:npx:noe);
for i=1:npx
for j=1:npy
k=(j-1)*npx+i;
trov=0;
k1=0;
while k1<=noe & trov==0
k1=k1+1;
if(ctx(i)==xyl(k1,1) & cty(j)==xyl(k1,2))
p(k)=k1;trov=1;
end
end
end
end

