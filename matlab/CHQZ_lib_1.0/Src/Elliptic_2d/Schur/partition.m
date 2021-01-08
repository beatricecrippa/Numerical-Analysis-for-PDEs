function [D]=partition(Rgamma)
% PARTITION Computes the diagonal weighting matrix D relative to interface unknowns
%
% [D]=partition(Rgamma)
%
%  Computes the diagonal weighting matrix D relative to the interface
%  unknowns: 
%  D_ii=1/n_i   where n_i is the number of subdomains x_i belongs to.
%  (see CHQZ3, pag. 397)
%
% Input : Rgamma = matrix contructed in cosrgam.m
%
% Output : D, column array of size ngamma (global number of interface nodes)
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
D=zeros(ngamma,1);
for i=1:ngamma
D(i)=max(Rgamma(:,i));
end
