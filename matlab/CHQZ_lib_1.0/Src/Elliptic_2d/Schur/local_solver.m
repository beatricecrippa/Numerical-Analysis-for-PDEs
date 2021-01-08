function[un]=local_solver(Amm,AGm,Lmm,LGG,novi,nvli,nov,novg,lint,lgamma,ugamma,f,ub);
% LOCAL_SOLVER  Solution of local problems after knowledge of u on the interface
%
% [un]=local_solver(Amm,AGm,Lmm,LGG,novi,nvli,nov,novg,lint,lgamma,ugamma,f,ub);
%
% Input: Amm = structure of local matrices (A_{m,m})^{-1}
%        AGm = structure of local matrices A_{\Gamma_m,m}
%        Lmm = structure containing maps from internal{\Omega_m to \Omega_m
%                set in schur_local.m
%        LGG = structure with maps from \Gamma_i to \Gamma
%                set in schur_local.m
%        novi = 2-indexes array of size (max(nvli),ne), computed in cosnovi
%        nvli = column array. nvli(ie) is the number of nodes of \Omega_ie
%        internal to Omega.
%        nov = 2-index array of local to global map, 
%                size(nov)=[max(npdx*npdy),ne]
%        novg = restriction map from (\cup_i \Omega_i) to \Gamma
%                set in cosnovg.m
%        lint =  list of those nodes of Omega (close) which
%                are internal to Omega.
%        lgamma = list of internal nodes which are on the interface between
%                spectral elements
%        ugamma = u on the interface  Gamma
%        f = r.h.s.
%        ub = solution on \partial\Omega
%
% Output: un = solution in Omega
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


ne=length(nvli);
noei=length(lint);
un_i=zeros(noei,1);
for ie=1:ne
    nn=novi(1:nvli(ie),ie);
    un_loc=zeros(nvli(ie),1);
    f_loc=f(nn);
    f_loc_i=f_loc(Lmm{ie})-(AGm{ie})'*ugamma(novg(LGG{ie},ie));
un_loc(Lmm{ie})=Amm{ie}*f_loc_i;
un_i(nn)=un_i(nn)+un_loc;
end
un_i(lgamma)=ugamma;
un=ub;
un(lint)=un_i;
return
