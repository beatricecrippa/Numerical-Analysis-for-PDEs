function [novi,nvli]=cosnovi(nov,ifro,lint)
% COSNOVI : constructs matrix novi realizing operator R_m
%
%   R_m is the restriction operator from the vector 
%   of coefficient unknowns related to the nodes of Omega (closed) to the
%   vector of coefficient unknowns related to the nodes of Omega_m (closed)
%
% Input: nov = 2-index array of local to global map, 
%                size(nov)=[max(npdx*npdy),ne] 
%        ifro = column array of length noe=nov(npdx*npdy,ne): 
%            if (x_i,y_i) is internal to Omega then ifro(i)=0,
%            if (x_i,y_i) is on \partial\Omega then ifro(i)=1,
%        lint =  list of those nodes of Omega (close) which
%                are internal to Omega.
%
% Output: novi = 2-indexes array of size (max(nvli),ne), computed in cosnovi
%         nvli = column array. nvli(ie) is the number of nodes of \Omega_ie
%                internal to Omega.
%

%   Written by Paola Gervasio
%   $Date: 2007/04/01$




[mn,ne]=size(nov);
nvli=zeros(ne,1);
nint=length(lint);
novi=zeros(mn,ne); 

for ie=1:ne
k=0;
for i=1:mn
ip=nov(i,ie);
if(ifro(ip)~=1)
k=k+1; trov=0;
j=1;
while j<=nint & trov==0
if(ip==lint(j))
novi(k,ie)=j; trov=1;
end
j=j+1;
end
end
end
nvli(ie)=k;
end

return
