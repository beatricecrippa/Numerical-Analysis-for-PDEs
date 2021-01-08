function [dlnp1]=derpol_legendre(n,ln,dln,dlnm1)
% DERPOL_LEGENDRE  recursive construction of Legendre polynomials 1st derivative
%            formula (2.3.19), pag. 77, CHQZ2
%
%   [dlnp1]=derpol_legendre(n,ln,dln,dlnm1)
%
% Input: n = polynomial degree
%        ln = character expression of L_n(x)
%        dln = character expression of (L_n)'(x)
%        dlnm1 = character expression of (L_{n-1})'(x)
%
% Output: dlnp1 = character expression of (L_{n+1})'(x)
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$


ns=num2str(n);
dlnp1=['(2*',ns,'+1)/(',ns,'+1).*(',ln,'+x.*(',dln,'))-',ns,'/(',ns,'+1).*(',dlnm1,')'];
return

