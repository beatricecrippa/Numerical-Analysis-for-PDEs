% FIG4_10  Script to produce  Fig 4.10, pag. 203 CHQZ2
% 
% Legendre collocation first-derivative eigenvalues computation and plot
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$


N=[8,16,32,64];
    nu=1;
    pbl=0;
% left
for nx=N
    [d,A]=lgl_eig(nx,nu,pbl);
fig=figure(...,
    'Name',['Fig. 4.10 N=',num2str(nx)],...
    'Visible','on');
plot(real(d),imag(d),'ko'); grid on
xlabel('Real'); ylabel('Imag');
axis equal
legend(['N=',num2str(nx)]);
end
