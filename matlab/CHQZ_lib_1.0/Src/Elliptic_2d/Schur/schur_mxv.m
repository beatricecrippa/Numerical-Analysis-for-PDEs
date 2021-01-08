function [v]=schur_mxv(x,AGG,Amm,AGm,LGG,novg,ne);
% SCHUR_MXV Computes matrix vector product, where the matrix is the Schur compl.
%
% Computes matrix vector product v=A*x, where matrix A is the Schur complement 
% matrix
%
%   [v]=schur_mxv(x,AGG,Amm,AGm,LGG,novg,ne);
%
% Input : x = column array, to perform v=A*x
%         AGG = structure of local matrices A_{\Gamma_m,\Gamma_m}
%         Amm = structure of local matrices (A_{m,m})^{-1}
%         AGm = structure of local matrices A_{\Gamma_m,m}
%         LGG = structure with maps from \Gamma_i to \Gamma
%                 set in schur_local.m
%         novg = restriction map from (\cup_i \Omega_i) to \Gamma
%                 set in cosnovg.m
%         ne = number of spectral elements
%
% Output : v= column array v=A*x, where A is the Schur complement matrix
%
% References: CHQZ3 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Evolution to Complex Geometries
%                     and Applications to Fluid DynamicsSpectral Methods"
%                    Springer Verlag, Berlin Heidelberg New York, 2007.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$


 
n=length(x);
v=zeros(n,1);

for ie=1:ne
lg=LGG{ie};
ngammam=length(lg);
agg=AGG{ie};agm=AGm{ie};amm=Amm{ie};
xloc=x(novg(lg,ie));
vloc=agg*xloc-agm*(amm*(agm'*xloc));
v(novg(lg,ie))=v(novg(lg,ie))+vloc;
end



