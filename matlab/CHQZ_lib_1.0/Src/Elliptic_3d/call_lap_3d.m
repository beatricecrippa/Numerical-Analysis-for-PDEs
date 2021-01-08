% CALL_LAP_3D Script for pre- and post-processing lap_3d
%

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

%
[uex,uex_x,uex_y,uex_z,ff,gam]=setfun_lap_3d;  % -Delta u + gam *u =f
%
xa=-1;xb=1;   % Omega=(xa,xb) x (ya,yb) x (za,zb)
ya=-1;yb=1;
za=-1;zb=1;
errinf=[];errh1=[];errl2=[];gammax=[]; gammay=[]; gammaz=[];
H=[];
for nex=4;
    ney=nex;  % decomposition of Omega in nex * ney * nez parallelepipeds
    nez=nex; 
    H=[H;2/nex];
for nx=2:4      % polynomial degree in each element along x-direction
    ny=nx;    % polynomial degree in each element along y-direction
    nz=nx;    % polynomial degree in each element along y-direction
    param=zeros(20,1);  
    param(1)=1;    % 1=SEM-NI,  
    param(2)=0;    % 0=no reordering, 1=CM ordering, 2=AMD ordering
    param(3)=2;    % 1= solve linear system by \
                   % 2= compute extrema eigenvalues of A
                   % 3 solve by Schur complement
                   % 4= compute extrema eigenvalues of the Schur complement
                   % 5= plot eigenvalues of matrix A
                   % 6: Compute extrema eigenvalues of  M^{-1}A
    param(4)=1;    % 1 computes errors, 0 no
    param(5)=0;    % 0 = errors with exact norms, 1= errors with discrete norms
    param(6)=nx*2; % n. of nodes for LG quadrature formulas
    param(7)=0;    % 0 =absolute errors, 1=relative errors
    param(8)=0;    % 0 no plot, 1 mesh, 2 surf, 3 contour
    param(9)=(nx+1); % nodes used to plot numerical solution
%     gammax=[0.2]; gammay=[0.2]; % if SEM decomposition is not uniform: 
%     gammaz=[0.2];            % they are the arrays with intefaces positions

% call lap_3d

[xy,un,D,param]=lap_3d(xa,xb,ya,yb,za,zb,gam,...
          uex,uex_x,uex_y,uex_z,ff,nex,nx,ney,ny,nez,nz,gammax,gammay,gammaz,...
          param);

% post

if (param(3)==1 | param(3)==3)
fprintf('nx=%d,nex=%d,err_inf=%11.4e, err_h1=%11.4e,err_l2=%11.4e \n',...
    nx,nex,param(29),param(30),param(31))
errinf=[errinf;param(29)];
errh1=[errh1;param(30)];
errl2=[errl2;param(31)];

elseif(param(3)==2)
fprintf('nx=%d,nex=%d,lam_max(A)=%11.4e, lam_min(A)=%11.4e, k(A)=%11.4e \n',...
    nx,nex,param(23),param(24),abs(param(23))/abs(param(24)));
elseif(param(3)==4)
fprintf('nx=%d,nex=%d,lam_max(S)=%11.4e, lam_min(S)=%11.4e, k(S)=%11.4e \n',...
    nx,nex,param(25),param(26),abs(param(25))/abs(param(26)));
elseif(param(3)==6)
fprintf('nx=%d,nex=%d,lam_max(M^{-1}A)=%11.4e, lam_min(M^{-1}A)=%11.4e, k(M^{-1}A)=%11.4e \n',...
    nx,nex,param(27),param(28),abs(param(27))/abs(param(28)));
end

end
end
