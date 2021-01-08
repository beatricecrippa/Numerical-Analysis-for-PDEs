function [etak]=pol_modal(k,lk,lkm2)
% POL_MODAL  recursive construction of modal basis function 
%            formula (2.3.31), pag. 82, CHQZ2
%   [etak]=pol_modal(k,lk,lkm2)
%            called by plot_modal
%
% Input: k = polynomial degree
%        lk = character expression of L_k(x)
%        lkm2 = character expression of L_{k-2}(x)
%
	% Output: etak = character expression of \eta_k(x)
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$


ks=num2str(k);
etak=['-1/sqrt(2*(2*',ks,'-1))*((',lk,')-(',lkm2,'))'];
return

