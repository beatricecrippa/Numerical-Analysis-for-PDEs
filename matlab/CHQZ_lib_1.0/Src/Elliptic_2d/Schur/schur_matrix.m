function [Sigma,PNN]=schur_matrix(ifro,nov,wx,dx,jacx,...
    wy,dy,jacy,nvli,gam,novg,lint,lgamma,D,Rgamma,param);
% SCHUR_MATRIX Computes  Schur complement matrix Sigma and its preconditioner
%
% Computes  Schur complement matrix Sigma and, if param(1)~=1, 
% the preconditioner for Sigma
%
%    [Sigma,PNN]=schur_matrix(ifro,nov,wx,dx,jacx,...
%    wy,dy,jacy,nvli,gam,novg,lint,lgamma,D,Rgamma,param);
%
% Input:
%         ifro = column array of length noe=nov(npdx*npdy,ne):
%            if (x_i,y_i) is internal to Omega then ifro(i)=0,
%            if (x_i,y_i) is on \partial\Omega then ifro(i)=1,
%        nov = local -global map, previously generated by cosnov_2d
%        wx = npdx LGL weigths in [-1,1], previously generated by xwlgl
%        dx =first derivative LGL matrix (by calling dx=derlgl(x,npdx))
%        jacx = array (length(jacx)=ne); jacx(ie)= (x_V2_ie-x_V1_ie)/2
%        wy = npdy LGL weigths in [-1,1],
%            (produced by calling [y,wy]=xwlgl(npdy))
%        dy =first derivative LGL matrix (by calling dy=derlgl(y,npdy))
%        jacy = array (length(jacy)=ne); jacy(ie)= (y_V3_ie-y_V1_ie)/2
%        nvli = column array with number of internal nodes of each
%               subdomain
%        novg = restriction map from (\cup_i \Omega_i) to \Gamma
%                 set in cosnovg.m
%        lint= list of internal nodes
%        lgamma= list of internal nodes which are on the interface between
%                spectral elements
%        D, column array of size ngamma (global number of interface nodes)
%        Rgamma = matrix contructed in cosrgam.m
%        param = array of parameters
%        
% Output: Sigma = Schur complement matrix
%         PNN = empty array if param(1)==1,
%               Neumann-Neumann preconditioner for Sigma, if param(1)==2,
%               bal Neumann-Neumann preconditioner for Sigma, if param(1)==3
%       
%
% References: CHQZ3 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Evolution to Complex Geometries 
%                     and Applications to Fluid DynamicsSpectral Methods"
%                    Springer Verlag, Berlin Heidelberg New York, 2007.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

[ldnov,ne]=size(nov);
noe=nov(ldnov,ne); 
npdx=length(wx); npdy=length(wy); mn=npdx*npdy;
n=length(lgamma);
Sigma=sparse(n,n);PSH=[];
if(param(1)==1)
    PNN=[];
else
    PNN=sparse(n,n);
end

[wx1,wy1]=meshgrid(wx,wy); w_ie=wx1.*wy1; w_ie=w_ie'; w_ie=w_ie(:); clear wx1 wy1;
if(param(1)==1) % P=I 
    
for ie=1:ne
Al=sparse(mn,mn);
[Al]=stiff_2d_sp(wx,dx,jacx(ie),wy,dy,jacy(ie));
jac=jacx(ie)*jacy(ie);

Al=Al+jac*gam*spdiags(w_ie,0,mn,mn);
ifro_l=ifro(nov(1:mn,ie));
[lbor,lint]=liste1(ifro_l);

Ali=Al(lint,lint);

ifroi=ifro_l(lint);
[lbor,lint,lintint,lg]=liste1(ifroi);


Amm=Ali(lintint,lintint); 
Agg=Ali(lg,lg);
Agm=Ali(lg,lintint);
Sigma_m=Agg-Agm*inv(Amm)*Agm';
Sigma(novg(lg,ie),novg(lg,ie))=Sigma(novg(lg,ie),novg(lg,ie))+Sigma_m;
end

elseif(param(1)==2)  %P= Neumann -Neumann preconditioner

for ie=1:ne
[Al]=stiff_2d_sp(wx,dx,jacx(ie),wy,dy,jacy(ie));
jac=jacx(ie)*jacy(ie);

Al=Al+jac*gam*spdiags(w_ie,0,mn,mn);

ifro_l=ifro(nov(1:mn,ie));
[lbor,lint]=liste1(ifro_l);
Ali=Al(lint,lint);
ifroi=ifro_l(lint);
[lbor,lint,lintint,lg]=liste1(ifroi);

Amm=Ali(lintint,lintint); 
Agg=Ali(lg,lg);
Agm=Ali(lg,lintint);

if nvli(ie)==mn
% Assembling preconditioner on the element ie.
% Mass/H^2 is added if the boundary of element ie  and boundary of 
% Omega have null intersection
% ww is not multiplied by jac, since Mass/H^2= diag(ww*jac)/jac
 Amp=Ali+spdiags(w_ie,0,mn,mn);
else
Amp=Ali;
end

Sigma_m=Agg-Agm*inv(Amm)*Agm';
Sigma(novg(lg,ie),novg(lg,ie))=Sigma(novg(lg,ie),novg(lg,ie))+Sigma_m;
clear Sigma_m;

Sigmap_m=inv(Amp(lg,lg)-Amp(lg,lintint)*inv(Amp(lintint,lintint))*Amp(lintint,lg));
D_m=D(novg(lg,ie));
PNN(novg(lg,ie),novg(lg,ie))=...
    PNN(novg(lg,ie),novg(lg,ie))+diag(D_m)*Sigmap_m*diag(D_m);
end

elseif(param(1)==3)  %P= balancing Neumann -Neumann preconditioner

for ie=1:ne
[Al]=stiff_2d_sp(wx,dx,jacx(ie),wy,dy,jacy(ie));
jac=jacx(ie)*jacy(ie);

Al=Al+jac*gam*spdiags(w_ie,0,mn,mn);

ifro_l=ifro(nov(1:mn,ie));
[lbor,lint]=liste1(ifro_l);
Ali=Al(lint,lint);
ifroi=ifro_l(lint);
[lbor,lint,lintint,lg]=liste1(ifroi);

Amm=Ali(lintint,lintint); 
Agg=Ali(lg,lg);
Agm=Ali(lg,lintint);

if nvli(ie)==mn
% Assembling preconditioner on the element ie.
% Mass/H^2 is added if the boundary of element ie  and boundary of 
% Omega have null intersection
% ww is not multiplied by jac, since Mass/H^2= diag(ww*jac)/jac
 Amp=Ali+spdiags(w_ie,0,mn,mn);
else
Amp=Ali;
end
nl=novg(lg,ie);
Sigma_m=Agg-Agm*inv(Amm)*Agm';
Sigma(nl,nl)=Sigma(nl,nl)+Sigma_m;
clear Sigma_m;

Sigmap_m=inv(Amp(lg,lg)-Amp(lg,lintint)*inv(Amp(lintint,lintint))*Amp(lintint,lg));
D_m=D(nl);
PNN(nl,nl)=PNN(nl,nl)+diag(D_m)*Sigmap_m*diag(D_m);
end

% A_H=Rgamma*Sigma*Rgamma'
% Sigma_H^+

PSH=Rgamma'*pinv(full(Rgamma*Sigma*Rgamma'))*Rgamma;
PNN=(speye(n)-PSH*Sigma)*PNN*(speye(n)-Sigma*PSH)+PSH;

end  % end loop on ie
return

