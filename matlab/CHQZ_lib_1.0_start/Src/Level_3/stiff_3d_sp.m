function [A]=stiff_3d_sp(wx,dx,jacx,wy,dy,jacy,wz,dz,jacz);
% STIFF_3D_SP Computes 3D local stiffness SEM matrix associated to (nabla(phi_j), nabla(phi_i))_N
%
%     (nabla(phi_j), nabla(phi_i))_N 
%
%
%    [A]=stiff_3d_sp(wx,dx,jacx,wy,dy,jacy,wz,dz,jacz);
%        produces the matrix
%        A_{ij}=(nabla(phi_j), nabla(phi_i))_N 
%        of size
%        (mn,mn) where mn=npdx*npdy is the local number of d.o.f.
%
% Input : 
%         wx = npdx LGL weigths in [-1,1],
%            (produced by calling [x,wx]=xwlgl(npdx))
%         dx =first derivative LGL matrix (by calling dx=derlgl(x,npdx))
%         jacx = array (length(jacx)=ne); jacx(ie)= (x_V2_ie-x_V1_ie)/2
%         wy = npdy LGL weigths in [-1,1],
%            (produced by calling [y,wy]=xwlgl(npdy))
%         dy =first derivative LGL matrix (by calling dy=derlgl(y,npdy))
%         jacy = array (length(jacy)=ne); jacy(ie)= (y_V3_ie-y_V1_ie)/2
%         wz = npdz LGL weigths in [-1,1],
%            (produced by calling [z,wz]=xwlgl(npdz))
%         dz =first derivative LGL matrix (by calling dzderlgl(z,npdz))
%         jacz= array (length(jacy)=ne); jacz(ie)= (y_V5_ie-y_V1_ie)/2
%         ww = column array with local weigths, length: npdx*npdy*npdz
%
% Output: A = matrix (npdx*npdy*npdz,npdx*npdy*npdz)
%
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$


npdx=length(wx);
npdy=length(wy);
npdz=length(wz);
jac=jacx*jacy*jacz;
jacyx=jacy/jacx;
jacxy=jacx/jacy;
mn=npdx*npdy;
ldnov=mn*npdz;
A=sparse(ldnov,ldnov);

% * dx' * diag(wx) * dx
sx = dx' * diag(wx) * dx;

% * dy' * diag(wy) * dy
sy = dy' * diag(wy) * dy;

% * dz' * diag(wz) * dz
sz = dz' * diag(wz) * dz;

for li=1:npdz; 
    il=(li-1)*mn;
    
for ki=1:npdy; 
    ik=(ki-1)*npdx;
coefx=wy(ki)*wz(li)*jacz*jacy/jacx;

for i=1:npdx; 
ii=il+ik+i;
coefy=wx(i)*wz(li)*jacz*jacx/jacy;
coefz=wx(i)*wy(ki)*jacy*jacx/jacz;

i1=il+ik+1;
i2=il+ik+npdx;
A(ii,i1:i2)=A(ii,i1:i2)+sx(i,1:npdx)*coefx;

i1=il+i;i2=il+mn;
A(ii,i1:npdx:i2)=A(ii,i1:npdx:i2)+sy(ki,1:npdy)*coefy;

i1=ik+i;
A(ii,i1:mn:ldnov)=A(ii,i1:mn:ldnov)+sz(li,1:npdz)*coefz;
end
end
end


% A=sparse(ldnov,mn);A=0;
% B=sparse(mn,mn);B=0;
% for ki=1:npdy
%     inde=((ki-1)*npdx+1:ki*npdx);
%     A(inde,inde)=dx'*(diag(ww(inde))*dx)*jacyx;
% end
% 
% for i=1:npdx
%     inde=(i:npdx:mn);
%     B(inde,inde)=dy'*(diag(ww(inde))*dy)*jacxy;
% end
% A=A+B;
% 
