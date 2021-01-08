function [nov_p1,nov_p1_g]=nov_p1_1d(nov);
% NOV_P1_1D  P1-FEM Local to global map
%
% [nov_p1,nov_p1_g]=nov_p1_1d(nov);
%
% INPUT: nov = map mesh_PN in Omega_ie -----> global mesh in Omega
%
% OUTPUT:   nov_p1 = map: mesh_P1 -----> mesh P_N in Omega_ie
%           nov_p1_g = map: mesh_P1  -----> global mesh in Omega

[npdx,ne]=size(nov);
noe=nov(npdx,ne);
nx=npdx-1;
nov_p1=zeros(2,nx);
nov_p1_g=zeros(2,nx);
for ie=1:nx
    nov_p1(1:2,ie)=[ie,ie+1];
    nov_p1_g(1:2,ie)=nov([ie,ie+1]);
end
return
