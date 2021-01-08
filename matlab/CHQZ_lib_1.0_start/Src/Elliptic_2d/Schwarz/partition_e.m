function [p_unity]=partition_e(nove,nvle,noe)
%  PARTITION_E   Unity partition for the extended mesh
%
%  [p_unity]=partition_e(nove,nvle,noe)
%
% Input : nove = 2-index array of "extended local" to global map, 
%        nvle = 2-index array:
%                 nvle (:,1)=number of nodes of the extendend elements
%                 nvle (:,2)=number of Q1 elements of the extendend elements
%                            along x-direction
%                 nvle (:,3)=number of Q1 elements of the extendend elements
%                            along y-direction
%        noe = global number of nodes
%
% Output : p_unity = Unity partition for the extended mesh

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

[ldnov,ne]=size(nove);
p_unity=zeros(noe,1);
for ie=1:ne
for i=1:nvle(ie,1);
    ii=nove(i,ie);
    p_unity(ii)=p_unity(ii)+1;
end
end
p_unity=1./p_unity;
return
