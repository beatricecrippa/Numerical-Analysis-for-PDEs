addpath Level_0
addpath Level_2
addpath Elliptic_2d

% CALL_LAP_2D Script for pre- and post- processing lap_2d
% Problem's data for - Delta u + gamma u = f
 
uex  = @(x,y) sin(pi*x).*sin(pi*y);   %Analytical solution
uexx = @(x,y) pi*cos(pi*x).*sin(pi*y);   %Gradient Analytical solution (x)
uexy = @(x,y) pi*sin(pi*x).*cos(pi*y);   %Gradient Analytical solution (y)
ff   = @(x,y) 2*pi^2*sin(pi*x).*sin(pi*y);   %Forcing term
g    = @(x,y) 0*x.*y;   %Dirichlet condition
h    = @(x,y) 0*x.*u;   %Neumann condition
gam  = 0;                 %gamma value

%Omega = (xa,xb) x (ya,yb)
xa = 0; xb = 1;  
ya = 0; yb = 1;
 
% boundary conditions on the sides of \partial\Omega
% numeration is 1-bottom, 2-right, 3-up, 4-left
% n = neumann, d = dirichlet
cb = ['dddd'];     



%%

err_inf = [];
err_H1 = [];
err_L2 = [];
h = [];
N = [];

for nex = [2 4 8 16 18 20]                 %[2,4,8,10,12,14,16,18,20];
    ney = nex;               % decomposition of Omega in nex x ney rectangles
    
    h = [h, 1/nex];
    
    erriN = [];
    errHN = [];
    errLN = [];
    
    for nx = [1 2 3 4 5]               % polynomial degree in each element along x-direction
        ny = nx;             % polynomial degree in each element along y-direction
        param=zeros(20,1);
        param(1)=1;       % 1=SEM-NI,  2= Patching
        param(2)=0;       % 0=no reordering, 1=CM ordering, 2=AMD ordering
        param(3)=1;       % 1= solve linear system by Choleski fact.
                          % 2= compute extrema eigenvalues of A
                          % 3= solve by Schur complement
                          % 4= compute extrema eigenvalues of the Schur complement
        param(4)=1;       % computes errors
        param(5)=0;       % 0 exact norms, 1= discrete norms
        param(6)=nx*2;    % nq for LG quadrature formulas
        param(7)=0;       % 0 =absolute errors, 1=relative errors
        param(8)=0;       % 0 no plot, 1 mesh, 2 surf, 3 contour
        param(9)=(nx+1);  % nodes used to plot numerical solution
        gammax=[]; gammay=[]; % if SEM decomposition is not uniform:
        % they are the arrays with intefaces positions
        
        % call lap_2d
        
        [xy,un,D,param]=lap_2d(xa,xb,ya,yb,gam,uex,uexx,uexy,ff,g,h,cb,...
            nex,nx,ney,ny,gammax,gammay,param);
        
        % output
        if (param(3)==1 || param(3)==3)
            fprintf('nx=%d,nex=%d,err_inf=%11.4e, err_h1=%11.4e,err_l2=%11.4e \n',...
                nx,nex,param(29),param(30),param(31))
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

        N = [N, nx];
        
        erriN = [erriN; param(29)];
        errHN = [errHN; param(30)];
        errLN = [errLN; param(31)];
    

    end
    
    err_inf = [err_inf, erriN];
    err_H1 = [err_H1, errHN];
    err_L2 = [err_L2, errLN];
end

figure(1);
hs = subplot(2,1,1);
loglog(h,h.^(N(3)+1),'-+b','Linewidth',2); hold on; grid on; 
loglog(h,err_L2(3,:),'-or','Linewidth',2); hold on; 
legend(sprintf('h^%i',N(3)+1),'||.||_{L^2}');
title ('Errors'); ylabel('L^2-error'); xlabel('h');
hs.FontSize = 12;

hs = subplot(2,1,2);
loglog(h,h.^N(3),'-+b','Linewidth',2); hold on; grid on;
loglog(h,err_H1(3,:),'-or','Linewidth',2); hold on;
legend(sprintf('h^%i',N(3)),'||.||_{H^1}');
ylabel('H^1-error'); xlabel('h');
hs.FontSize = 12;

