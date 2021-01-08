function [z]=preconnl(r,ne,nvli,novg,D,LGG,Am);
% SCHUR_PRECONNL:  Solves the sistem (P^{NN})^{-1} z=r where P^{NN} is the Neumann-Neumann preconditioner for Schur compl.  
%
%   Neumann-Neumann preconditioner for Schur
% Complement    (Algorithm 6.4.3., pag. 398 CHQZ3)
%                   Solves the sistem (P^{NN})^{-1} z=r
%
% Input: r= r.h.s, column array of length ngamma
%        ne = number of subdomains
%        nvli = column array with number of internal nodes of each
%               subdomain
%        novg = restriction map from (\cup_i \Omega_i) to \Gamma
%                 set in cosnovg.m
%        D = diagonal weighting matrix (column array of lenght ngamam)
%        LGG = structure with maps from \Gamma_i to \Gamma
%                 set in schur_local.m
%        Am  = structure of choleski factors of local matrices 
%              \tilde A_m= A_mm+ M_m/H_m^2 (H_m is the size of the
%              subdomain of index m, while M_m si the local mass matrix)
%
% Ouput: z = solution, column array of length ngamma
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


ngamma=length(r);

z=zeros(ngamma,1);
for ie=1:ne
% 
lgamma=LGG{ie};
mn=nvli(ie);
% compute D_ie (restriction of global weighting matrix  D to Gamma_ie)
Die=D(novg(lgamma,ie));
f=zeros(mn,1);
f(lgamma)=Die.*r(novg(lgamma,ie));
R=Am{ie};
zloc=R\(R'\f);
z(novg(lgamma,ie))=z(novg(lgamma,ie))+Die.*zloc(lgamma);
end
return
