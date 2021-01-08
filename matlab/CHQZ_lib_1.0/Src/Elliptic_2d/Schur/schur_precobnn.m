function [z]=schur_precobnn(r,ne,nvli,novg,D,LGG,Am,Rgamma,PSH,AGG,Amm,AGm);
% SCHUR_PRECOBNN:   Solves the sistem (P_B^{NN})^{-1} z=r where P_B^{NN} is the Balancing Neumann-Neumann preconditioner for Schur compl.
%    (Algorithm 6.4.4., pag. 399 CHQZ3)
%                   Solves the sistem (P_B^{NN})^{-1} z=r
%
% Input: r= r.h.s, column array of length ngamma
%        ne = number of subdomains
%        nvli = column array with number of internal nodes of each
%        subdomain
%        novg = restriction map from (\cup_i \Omega_i) to \Gamma
%                 set in cosnovg.m
%        D = diagonal weighting matrix (column array of lenght ngamam)
%        LGG = structure with maps from \Gamma_i to \Gamma
%                 set in schur_local.m
%        Am  = structure of choleski factors of local matrices 
%              \tilde A_m= A_mm+ M_m/H_m^2 (H_m is the size of the
%              subdomain of index m, while M_m si the local mass matrix)
%        Rgamma = matrix of size (ne, ngamma) (defined in CHQZ3, pag. 399)
%                  (R_\Gamma)_ij: =  1/n_j     if  x_j \in \Gamma_i
%                                 =  0         otherwise
%        PSH  = pseudoinverse of Sigma_H, i.e. Rgamma^T A_H^+ Rgamma
%        AGG = structure of local matrices A_{\Gamma_m,\Gamma_m}
%        Amm = structure of local matrices (A_{m,m})^{-1}
%        AGm = structure of local matrices A_{\Gamma_m,m}
%
% Ouput: z = solution, column array of length ngamma
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


ra=PSH*r;

y=r-schur_mxv(ra,AGG,Amm,AGm,LGG,novg,ne);
z=schur_preconnl(y,ne,nvli,novg,D,LGG,Am);
y=schur_mxv(z,AGG,Amm,AGm,LGG,novg,ne);
z=z-PSH*y+ra;
return
