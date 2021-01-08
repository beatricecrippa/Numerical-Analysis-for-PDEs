function [Al,Ml,Mld]=matricesp1_1d(nu,beta,gam,jacx,param);
% MARICESP1_1D     P1 - Local mass and stiffness matrices on [-1,1]
%
%  [Al,Ml,Mld]=matricesp1_1d(jacx);
%
% Input: nu   = viscosity (constant>0)
%        beta  = coefficient of first order term (constant)
%        gam   = coefficient of zeroth order term (constant>=0)
%        jacx = jacobian of the map F:[-1,1] --- > (x_1,x_2)
%
% Output: Al  =local P1 matrix
%         Ml  = Local P1 mass matrix
%         Mld = local P1 mass matrix with numerical integration (trapezoidal
%        rule)
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$


% local P1 - mass matrix
Ml=[2 1; 1 2];
Ml=Ml/3*jacx;

% local P1 - stiffness matrix
Al=[1, -1; -1, 1];
Al=Al/(2*jacx);

% local P1 - first derivative matrix
Al1=[-1, 1; -1, 1];
Al1=Al1/2;

% local P1  - mass matrix with numerical integration 
Mld=spdiags(ones(2,1),0,2,2)*jacx;

if nu~=1
Al=nu*Al;
end
if beta~=0
    Al=Al+beta*Al1;
end
if gam~=0
    if param(1)==1 | param(1)==4 | param(1)==5
    Al=Al+gam*Ml;
    elseif param(1)==2 | param(1)==3
    Al=Al+gam*Mld;
    end
end
