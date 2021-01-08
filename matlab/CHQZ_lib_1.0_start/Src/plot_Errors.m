% plot convergence errors
% Ex 1 - h convergence
N = 3;
h = [0.5, 0.25, 0.125, 0.0625];
err_h1 = [1.7565e-03, 2.2164e-04, 2.7770e-05, 3.4733e-06];
err_l2 = [9.2104e-05, 5.8334e-06, 3.6579e-07, 2.2881e-08];

% ERROR PLOTS
figure(1);
hs = subplot(2,1,1);
loglog(h,h.^(N+1),'-+b','Linewidth',2); hold on; grid on; 
loglog(h,err_l2,'-or','Linewidth',2); hold on; 
legend(sprintf('h^%i',N+1),'||.||_{L^2}');
title ('Errors'); ylabel('L^2-error'); xlabel('h');
hs.FontSize = 12;

hs = subplot(2,1,2);
loglog(h,h.^N,'-+b','Linewidth',2); hold on; grid on;
loglog(h,err_h1,'-or','Linewidth',2); hold on;
legend(sprintf('h^%i',N),'||.||_{H^1}');
ylabel('H^1-error'); xlabel('h');
hs.FontSize = 12;

%%
% Ex 1 - N convergence
N = [2,3,4,5,6];
h = 0.25;
err_h1 = [1.0508e-02, 2.2164e-04, 3.4881e-06, 4.3806e-08, 4.5781e-10];
err_l2 = [4.0764e-04, 5.8334e-06, 7.0219e-08, 7.1558e-10, 6.2984e-12];

% ERROR PLOTS
figure(2);
hs = subplot(2,1,1);
semilogy(N,err_l2,'-or','Linewidth',2); hold on;  grid on;
legend('||.||_{L^2}');
title ('Errors'); ylabel('L^2-error'); xlabel('N');
hs.FontSize = 12;

hs = subplot(2,1,2);
loglog(N,err_h1,'-or','Linewidth',2); hold on; grid on;
legend('||.||_{H^1}');
ylabel('H^1-error'); xlabel('N');
hs.FontSize = 12;


%% Exercise 2
N = 1;
h = [0.5, 0.25, 0.125, 0.0625];
% err_h1 = [3.2606e+00, 2.1710e+00, 1.0589e+00, 5.2306e-01]; 
err_h1 = [7.5555e-01, 3.6957e-01, 1.8073e-01, 8.8866e-02];
figure(4); hs = gca;
loglog(h,h.^1,'-+b','Linewidth',2); hold on; grid on; 
loglog(h,err_h1,'-or','Linewidth',2); hold on; 
legend(sprintf('h^%i',1),'||.||_{H^1}');
title ('Errors N = 1'); ylabel('H^1-error'); xlabel('h');
hs.FontSize = 12;


N = 2; 
err_h1 = [9.3918e-02, 3.5882e-02, 1.4761e-02,6.2802e-03];
figure(5); hs = gca;
loglog(h,h.^1.5,'-+b','Linewidth',2); hold on; grid on; 
loglog(h,err_h1,'-or','Linewidth',2); hold on; 
legend(sprintf('h^{1.5}'),'||.||_{H^1}');
title ('Errors N = 2'); ylabel('H^1-error'); xlabel('h');
hs.FontSize = 12;


N = 3;
err_h1 = [3.3007e-02,1.4197e-02,6.1528e-03,2.6725e-03];
figure(6); hs = gca;
loglog(h,h.^1.5,'-+b','Linewidth',2); hold on; grid on; 
loglog(h,err_h1,'-or','Linewidth',2); hold on; 
legend(sprintf('h^{1.5}'),'||.||_{H^1}');
title ('Errors N = 3'); ylabel('H^1-error'); xlabel('h');
hs.FontSize = 12;


N = 4;
err_h1 = [1.7612e-02, 7.6341e-03,3.3157e-03,1.4416e-03];
figure(7); hs = gca;
loglog(h,h.^1.5,'-+b','Linewidth',2); hold on; grid on; 
loglog(h,err_h1,'-or','Linewidth',2); hold on; 
legend(sprintf('h^{1.5}'),'||.||_{H^1}');
title ('Errors N = 4'); ylabel('H^1-error'); xlabel('h');
hs.FontSize = 12;
