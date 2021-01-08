function [d,A]=lgl_eig(nx,nu,pbl);
% LGL_EIG   Computes eigenvalues of first/second order spectral derivative matrices: collocation/G-NI,  LGL nodes, 1D monodomain
%
%  [d,A]=lgl_eig(nx,nu,pbl);
%
% Input: nx = polynomial degree  
%        nu =  viscosity (constant > 0)
%        pbl = integer parameter to set the problem:
%              0 : compute eigenvalues of  A with
%              Au := du/dx    in (-1,1),  u(1)=0
%              strong form (collocation)
%
%            = 1 : compute eigenvalues of A with
%              Au := du/dx    in (-1,1),   u(1)=0    
%              weak form with integration by part
%
%            = 2 : compute generalized eigenvalues of A with
%              Au := du/dx    in (-1,1),   u(1)=0    
%              weak form with integration by part
%
%            = 30 : compute eigenvalues of A with
%              Au := -d^2u/dx^2    in (-1,1),   u(-1)=u(1)=0    
%              strong form
%
%            = 31 : compute eigenvalues of A with
%              Au := -d^2u/dx^2    in (-1,1),   u(-1)=u(1)=0    
%              weak form (G-NI)
% 
%            = 32 : compute generalized eigenvalues of A with
%              Au := -d^2u/dx^2    in (-1,1),   u(-1)=u(1)=0    
%              weak form (G-NI)
%
%            = 41 : compute eigenvalues of A with
%              Au := -d^2u/dx^2    in (-1,1),   u'(-1)=u'(1)=0    
%              weak form (G-NI)
% 
%            = 42 : compute generalized eigenvalues of A with
%              Au := -d^2u/dx^2    in (-1,1),   u'(-1)=u'(1)=0  
%              weak form (G-NI)
%
%            = 51 : compute eigenvalues of A with
%              Au := -nu d^2u/dx^2+du/dx    in (-1,1),   u(-1)=u(1)=0    
%              weak form (G-NI)
% 
%            = 52 : compute generalized eigenvalues of A with
%              Au := -nu d^2u/dx^2+du/dx    in (-1,1),   u(-1)=u(1)=0  
%              weak form (G-NI)
%
%            = 61 : compute eigenvalues of A with
%              Au := -nu d^2u/dx^2+du/dx    in (-1,1),   u(-1)=0, u'(1)=0    
%              weak form (G-NI)
% 
%            = 62 : compute generalized eigenvalues of A with
%              Au := -nu d^2u/dx^2+du/dx    in (-1,1),   u(-1)=0, u'(1)=0   
%              weak form (G-NI)
%
% Output: d = column array with eigenvalues of A
%         A = matrix required, of size (npdx,npdx)
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

%
npdx=nx+1;
[x,wx]=xwlgl(npdx); dx=derlgl(x,npdx);

% A

%              0 : compute eigenvalues of  A with
%              Au := du/dx    in (-1,1),  u(1)=0
%              strong form (collocation)
if pbl == 0
A=dx(1:nx,1:nx); 

elseif pbl == 1
%            = 1 : compute eigenvalues of A with
%              Au := du/dx    in (-1,1),   u(1)=0    
%              weak form with integration by part
A=-dx'*diag(wx);
A(1,1)=A(1,1)-1;

elseif pbl == 2
%            = 2 : compute generalized eigenvalues of A with
%              Au := du/dx    in (-1,1),   u(1)=0    
%              weak form with integration by part
A=-dx'*diag(wx);
A(1,1)=A(1,1)-1;
A=diag(1./wx)*A;

elseif pbl == 30
%            = 30 : compute eigenvalues of A with
%              Au := -d^2u/dx^2    in (-1,1),   u(-1)=u(1)=0    
%              strong form
A=-dx*dx; 
A=A(2:nx,2:nx);

elseif pbl == 31
%            = 31 : compute eigenvalues of A with
%              Au := -d^2u/dx^2    in (-1,1),   u(-1)=u(1)=0    
%              weak form (G-NI)
A=dx'*diag(wx)*dx;
A=A(2:nx,2:nx);

elseif pbl == 32
%            = 32 : compute generalized eigenvalues of A with
%              Au := -d^2u/dx^2    in (-1,1),   u(-1)=u(1)=0    
%              weak form (G-NI)
A=dx'*diag(wx)*dx;
A=diag(1./wx(2:nx))*A(2:nx,2:nx);
d=eig(A);

elseif pbl == 41
%            = 41 : compute eigenvalues of A with
%              Au := -d^2u/dx^2    in (-1,1),   u'(-1)=u'(1)=0    
%              weak form (G-NI)
A=dx'*diag(wx)*dx;
A=A(1:npdx,1:npdx);

elseif pbl == 42
%            = 42 : compute generalized eigenvalues of A with
%              Au := d^2u/dx^2    in (-1,1),   u'(-1)=u'(1)=0  
%              weak form (G-NI)
%
A=dx'*diag(wx)*dx;
A=diag(1./wx)*A;

elseif pbl == 51
%            = 51 : compute eigenvalues of A with
%              Au := -nu d^2u/dx^2+du/dx    in (-1,1),   u(-1)=u(1)=0    
%              weak form (G-NI)
A=dx'*diag(wx)*dx;
A=A(2:nx,2:nx);

elseif pbl == 52
%            = 52 : compute generalized eigenvalues of A with
%              Au := -nu d^2u/dx^2+du/dx    in (-1,1),   u(-1)=u(1)=0  
%              weak form (G-NI)
A=dx'*diag(wx)*dx;A=A*nu;
A=A-dx'*diag(wx);
A=diag(1./wx)*A;
A=A(2:nx,2:nx);

elseif pbl == 61
%            = 61 : compute eigenvalues of A with
%              Au := -nu d^2u/dx^2+du/dx    in (-1,1),   u(-1)=0, u'(1)=0    
%              weak form (G-NI)
A=dx'*diag(wx)*dx;A=A*nu;
A=A-dx'*diag(wx);
A(npdx,npdx)=A(npdx,npdx)+1;
A=A(2:npdx,2:npdx);

elseif pbl == 62
%            = 62 : compute generalized eigenvalues of A with
%              Au := -nu d^2u/dx^2+du/dx    in (-1,1),   u(-1)=0, u'(1)=0   
%              weak form (G-NI)
A=dx'*diag(wx)*dx;A=A*nu;
A=A-dx'*diag(wx);
A(npdx,npdx)=A(npdx,npdx)+1;
A=diag(1./wx)*A;
A=A(2:npdx,2:npdx);
end
d=eig(A);
