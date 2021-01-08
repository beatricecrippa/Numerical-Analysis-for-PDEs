% CALL_SCHWARZ_2D Script for pre and post processing schwarz_2d.
%

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

[uex,uex_x,uex_y,ff,g,h,gam]=setfun_lap_2d;  % -Delta u + gam *u =f
%
xa=-1;xb=1;   % Omega=(xa,xb) x (ya,yb)
ya=-1;yb=1;
cb='dddd'; % schwarz_2d works only if cb='dddd';
    param=zeros(20,1);  
    param(1)=2;    % 1:P=I, 2: P=P_as
    param(2)=2; % number of levels
    param(3)=2;   % 1=PCG, 2=PBicgstab
    param(4)=1;   % computes errors
fprintf('nx   nex   iter      res         err_inf          err_h1       err_l2\n')
for nex=4;
    ney=nex;  % decomposition of Omega in nex x ney rectangles
for nx=4:16  % polynomial degree in each element along x-direction
    ny=nx;     % polynomial degree in each element along y-direction
    param(5)=1;    % 0 exact norms, 1= discrete norms
    param(6)=nx*2;   % nq for LG quadrature formulas
    param(7)=1;    % 0 =absolute errors, 1=relative errors
    param(8)=0;    % 0 no plot, 1 mesh, 2 surf, 3 contour
    param(9)=(nx+1); % nodes used to plot numerical solution
    param(10)=1.d-12; % tolerance for pcg
    param(11)=400; % maxit for pcg
    gammax=[]; gammay=[]; % if SEM decomposition is not uniform: 
                          % they are the arrays with intefaces positions

   % call schur

[xy,un,param]=schwarz_2d(xa,xb,ya,yb,gam,...
          uex,uex_x,uex_y,ff,g,h,cb,nex,nx,ney,ny,gammax,gammay,param);
 
% output
fprintf('%d    %d     %d    %11.4e      %11.4e     %11.4e %11.4e \n',...
    nx,nex,param(21), param(22), param(25),param(26),param(27))

end
end
