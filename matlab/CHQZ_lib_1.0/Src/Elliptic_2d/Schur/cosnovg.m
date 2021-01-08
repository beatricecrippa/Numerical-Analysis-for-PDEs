function [novg]=cosnovg(xyi,noei,ifroi,lgamma,ldnov,novi,nvli);
% COSNOVG : constructs matrix novg realizing operator R_{\Gamma_m}
%
%   R_{\Gamma_m}: the restriction operator from the vector 
%   of coefficient unknowns related
%   to the nodes of Gamma to only those associated with Gamma_m=Gamma \cap
%   \partial\Omega_m.  (see CHQZ3, pag. 394)
%
% [novg]=cosnovg(xyi,noei,ifroi,lgamma,ldnov,novi,nvli);
%
% Input: xyi = 2-indexes array of coordinates of nodes internal to Omega
%        noei = number of nodes internal to Omega
%        ifroi =  restriction of ifro to nodes internal to Omega
%        lgamma =  list of those nodes of Omega (without boundary) which
%        belong to Gamma
%        ldnov, leading dimension of nov
%        novi = 2-indexes array of size (max(nvli),ne), computed in cosnovi
%        nvli = column array. nvli(ie) is the number of nodes of \Omega_ie
%        internal to Omega.
%
% Output: novg = 2-indexes array of size (ldnov,ne)
%         novg(i,ie)= 0 if node x_i of Omega_ie does not belong to Gamma
%                    = j if node x_i of Omega_ie is the node j of Gamma 
%
% References: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.
%             CHQZ3 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Evolution to Complex Geometries 
%                     and Applications to Fluid DynamicsSpectral Methods"
%                    Springer Verlag, Berlin Heidelberg New York, 2007.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$



ngamma=length(lgamma);
ne=length(nvli);
novg=zeros(ldnov,ne);
for ie=1:ne
mn=nvli(ie);
for i=1:mn
ig=novi(i,ie);
if(ifroi(ig)==-1)
x=xyi(ig,1); y=xyi(ig,2);
trov=0; j=0;
while j <=ngamma & trov==0
j=j+1; 
x1=xyi(lgamma(j),1); y1=xyi(lgamma(j),2);
if(x==x1 & y==y1)
trov=1; novg(i,ie)=j;
end
end
end
end
end

return
