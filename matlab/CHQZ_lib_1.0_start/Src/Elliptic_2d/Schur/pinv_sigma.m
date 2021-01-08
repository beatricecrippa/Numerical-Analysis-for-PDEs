function [PSH]=pinv_sigma(AGG,Amm,AGm,LGG,novg,Rgamma);
% PINV_SIGMA Computes the pseudoinverse of Sigma_H
%
% Computes the pseudoinverse of Sigma_H (Algorithm 6.4.4, pag. 399, CHQZ3) 
%
% Input : AGG = structure of local matrices A_{\Gamma_m,\Gamma_m}
%         Amm = structure of local matrices (A_{m,m})^{-1}
%         AGm = structure of local matrices A_{\Gamma_m,m}
%         LGG = structure with maps from \Gamma_i to \Gamma
%                 set in schur_local.m
%         novg = restriction map from (\cup_i \Omega_i) to \Gamma
%                 set in cosnovg.m
%         Rgamma = matrix of size (ne, ngamma) (defined in CHQZ3, pag. 399)
%                  (R_\Gamma)_ij: =  1/n_j     if  x_j \in \Gamma_i
%                                 =  0         otherwise
%
% Output: PSH  = pseudoinverse of Sigma_H, i.e. Rgamma^T A_H^+ Rgamma
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

[ne,ngamma]=size(Rgamma);
Sigma=zeros(ngamma,ngamma);
for ie=1:ne
lg=LGG{ie};
ngammam=length(lg); nl=novg(lg,ie);
agg=AGG{ie};agm=AGm{ie};amm=Amm{ie};
Sigma_loc=agg-agm*amm*agm';
Sigma(nl,nl)=Sigma(nl,nl)+Sigma_loc;
end

% A_H=Rgamma*Sigma*Rgamma'

PSH=pinv(full(Rgamma*Sigma*Rgamma'));
clear Sigma;
% Sigma_H^+

PSH=Rgamma'*PSH*Rgamma;
