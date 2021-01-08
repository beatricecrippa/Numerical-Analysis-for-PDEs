function [lnp1]=pol_legendre(n,ln,lnm1)
% POL_LEGENDRE  recursive construction of Legendre basis function
%            formula (2.3.19), pag. 77, CHQZ2
%   [lnp1]=pol_legendre(n,ln,lnm1)
%            called by plot_legendre
%
% Input: n = polynomial degree
%        ln = character expression of L_n(x)
%        lnm1 = character expression of L_{n-1}(x)
%
% Output: lnp1 = character expression of L_{n+1}(x)
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

ns=num2str(n);
lnp1=['(2*',ns,'+1)/(',ns,'+1).*x.*(',ln,')-',ns,'/(',ns,'+1).*(',lnm1,')'];
return

