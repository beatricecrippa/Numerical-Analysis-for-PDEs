function [Sigma]=schur_assemb(AGG,Amm,AGm,LGG,novg,...
    lint,lgamma,param);
% SCHUR_ASSEMB Assembles global Schur complement matrix
%
% Input: AGG =structure of local matrices A_{\Gamma_m,\Gamma_m}
%        Amm = structure of local matrices (A_{m,m})^{-1}
%        AGm = structure of local matrices A_{\Gamma_m,m}
%        Lmm = structure containing maps from internal{\Omega_m to \Omega_m
%                set in schur_local.m
%        LGG = structure with maps from \Gamma_i to \Gamma
%                set in schur_local.m
%        novi = 2-indexes array of size (max(nvli),ne), computed in cosnovi
%        nvli = column array. nvli(ie) is the number of nodes of \Omega_ie
%        internal to Omega.
%        nov = 2-index array of local to global map, 
%                size(nov)=[max(npdx*npdy),ne]
%        novg = restriction map from (\cup_i \Omega_i) to \Gamma
%                set in cosnovg.m
%        lint =  list of those nodes of Omega (close) which
%                are internal to Omega.
%        lgamma = list of internal nodes which are on the interface between
%                spectral elements
%        ugamma = u on the interface  Gamma
%        f = r.h.s.
%        ub = solution on \partial\Omega
%
% Output: un = solution in Omega
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


n=length(lgamma);
[ldnov,ne]=size(novg);
Sigma=zeros(n,n);

    % Schur complement 
for ie=1:ne
lg=LGG{ie};
ngammam=length(lg);
agg=AGG{ie};agm=AGm{ie};amm=Amm{ie};
Sigma_ie=agg-agm*amm*agm';
Sigma(novg(lg,ie),novg(lg,ie))=Sigma(novg(lg,ie),novg(lg,ie))+Sigma_ie;
end

if(param(1)==2)
    % (Neumann Neumann preconditioner)
for ie=1:ne
lg=LGG{ie};
ngammam=length(lg);
agg=AGG{ie};agm=AGm{ie};amm=Amm{ie};
Sigma_ie=agg-agm*amm*agm';
Sigma(novg(lg,ie),novg(lg,ie))=Sigma(novg(lg,ie),novg(lg,ie))+Sigma_ie;
end
end

return
