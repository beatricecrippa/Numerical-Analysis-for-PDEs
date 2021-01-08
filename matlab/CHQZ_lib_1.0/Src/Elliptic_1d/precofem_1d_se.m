function [AFE,MFE,MFEd]=precofem_1d_se(AFE,MFE,MFEd,ne,xy,nov,nu,beta,gam,param);
% PRECOFEM_1D_SE   P1 matrices (stiffness, mass, discrete mass) on macro spectral elements
%
% [AFE,MFE,MFEd]=precofem_1d_se(AFE,MFE,MFEd,ne,xy,nov,nu,beta,gam,param);
%
% Input:  AFE = stiffness P1 matrix
%         MFE = mass P1 matrix
%         MFEd = mass P1 matrix with numerical integration
%         ne =  number of P1 element inside spectral element Omega_ie
%         xy = mesh on the spectral element Omega_ie
%         nov = map mesh Q_n in Omega_ie ---> global mesh in Omega 
%         nu   = viscosity (constant>0)
%         beta  = coefficient of first order term (constant)
%         gam   = coefficient of zeroth order term (constant>=0)
%
% Output: AFE = stiffness P1 matrix 
%         MFE = mass P1 matrix
%         MFEd = mass P1 matrix with numerical integration
%
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

%    nov_q1 map: mesh_Q1 -----> mesh Q_N in Omega_ie
%    nov_q1_g map: mesh_Q1 -----> global mesh in Omega
%
    [nov_p1,nov_p1_g]=nov_p1_1d(nov);
    
for ie_p1=1:ne
    jacx=5.d-1*(xy(nov_p1(2,ie_p1))-xy(nov_p1(1,ie_p1)));
    [Al,Ml,Mld]=matricesp1_1d(nu,beta,gam,jacx,param);
    AFE(nov_p1_g(1:2,ie_p1),nov_p1_g(1:2,ie_p1))=AFE(nov_p1_g(1:2,ie_p1),nov_p1_g(1:2,ie_p1))+Al;
    MFE(nov_p1_g(1:2,ie_p1),nov_p1_g(1:2,ie_p1))=MFE(nov_p1_g(1:2,ie_p1),nov_p1_g(1:2,ie_p1))+Ml;
    MFEd(nov_p1_g(1:2,ie_p1),nov_p1_g(1:2,ie_p1))=MFEd(nov_p1_g(1:2,ie_p1),nov_p1_g(1:2,ie_p1))+Mld;
end
return
