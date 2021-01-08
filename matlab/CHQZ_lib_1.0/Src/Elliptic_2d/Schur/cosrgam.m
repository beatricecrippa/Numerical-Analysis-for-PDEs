function [Rgamma]=cosrgam(novg,LGG,ne,ngamma);
% COSRGAM  Computes  matrix R_Gamma for Schur complement preconditioners
%
%   [Rgamma]=cosrgam(novg,LGG,ne,ngamma);
%
% Computes  matrix R_Gamma of size (ne, ngamma), defined as follows:
%
% If \Gamma is the interface (uninion of interfaces between subdomains)
% and \Gamma_i :=\partial \Omega_i \cap \Gamma, then
% (R_\Gamma)_ij: =  1/n_j     if  x_j \in \Gamma_i
%                =  0         otherwise
%
%   where n_j = is the number of subdomains x_j belongs to
%
% (see CHQZ3, pag. 399)
%
% Matrix R_gamma is needed to construct the coarse operator A_H for 
% Balancing Neumann Neumann preconditioner (Schur complement matrix)
% and also to set the unity partition matrix D used in schur_preconn 
% and in preco_bnn
%
% Input : novg = restriction map from (\cup_i \Omega_i) to \Gamma
%                 set in cosnovg.m
%         LGG   = structure with maps from \Gamma_i to \Gamma
%                 set in schur_local.m
%         ne = number of subdomains
%         ngamma = number of nodes on Gamma
%
% Output : Rgamma
%
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


Rgamma=zeros(ne,ngamma);

for ie=1:ne
lgamma=LGG{ie};
Rgamma(ie,novg(lgamma,ie))=Rgamma(ie,novg(lgamma,ie))+1;
end

for i=1:ngamma
Rgamma(:,i)=Rgamma(:,i)/sum(Rgamma(:,i));
end

return
