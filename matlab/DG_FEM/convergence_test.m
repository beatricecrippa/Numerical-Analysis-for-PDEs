function [table_error,rates]= convergence_test(TestName,nRef)
%% [table_error,rates] = convergence_test(TestName,nRef)
%==========================================================================
% Error analysis varying the mesh size h 
%==========================================================================
% Example of usage: [table_error,rates] = convergence_test('Test1',[2:5])
%
%    INPUT:
%          test_name    : (string)  test case name, see dati.m
%          nRef         : (array int) vector for successive mesh refinement 
%
%    OUTPUT:
%          table_error  : (struct) containing the computed L^2, H^1 and DG errors
%          rates        : (struct) containing the computed rates


k = 1;
for i = 1:nRef
    [errors,solutions,femregion,Data]= main2D(TestName,i);
    E_L2(k) = errors.E_L2;
    E_H1(k) = errors.E_H1;
    E_DG(k) = errors.E_DG;
    ne(k) =femregion.ne;
    hh(k) = femregion.h;
    condNumA(k) = Data.condA;
    fprintf('End test %d\n',k);
    k = k + 1;
end

close all;

table_error = [ne',hh',E_L2',E_H1'];
rates = [log10(E_L2(2:end)./E_L2(1:end-1))', log10(E_H1(2:end)./E_H1(1:end-1))', log10(E_DG(2:end)./E_DG(1:end-1))']./log10(hh(2:end)/hh(1:end-1));

figure;
% hs = subplot(2,1,1);
p = femregion.degree;
loglog(hh,hh,'-+b','Linewidth',2);
hold on
loglog(hh,E_L2,'-or','Linewidth',2);
hold on
legend(sprintf('h'),'||.||_{L^2}');
title ('Errors');
ylabel('L^2-error')
xlabel('h');
hs.FontSize = 12;
% 
% hs = subplot(2,1,2);
% loglog(hh,hh.^(p),'-+b','Linewidth',2);
% hold on
% loglog(hh,E_H1,'-or','Linewidth',2);
% hold on
% loglog(hh,E_DG,'-pm','LineWidth',2);
% legend(sprintf('h^%i',p),'||.||_{H^1}','||.||_{DG}');
% ylabel('H^1-error')
% xlabel('h');
% hold off
% hs.FontSize = 12;
%  
%    
figure;
loglog(hh,1./hh.^2,'-+b','Linewidth',2);
xlabel('h');
ylabel('condition number');
hold on
loglog(hh,condNumA,'-or','Linewidth',2);
ha = gca; ha.FontSize = 12;
legend(sprintf('h^%d',2),'cond(A)');




   

 
   
   
   


