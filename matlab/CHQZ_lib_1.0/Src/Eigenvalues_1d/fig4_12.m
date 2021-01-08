% FIG4_12  Script to produce Fig 4.12 (top-left) and ,Fig 4.13 (top-left) pag. 204 CHQZ2
% 
% Legendre collocation/ G-NI / generalized G-NI first-derivative eigenvalues
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$



    nu=1;
fig=figure(...,
    'Name',['Fig. 4.12'],...
    'Visible','on');

nx=16 % or nx= ....
    pbl=0;
    [d,A]=lgl_eig(nx,nu,pbl);
plot(real(d),imag(d),'ko'); 
hold on
    pbl=1;
    [d,A]=lgl_eig(nx,nu,pbl);
plot(real(d),imag(d),'k*');
    pbl=2;
    [d,A]=lgl_eig(nx,nu,pbl);
plot(real(d),imag(d),'k>');
legend('LC','LG-NI','Gen LG-NI',2);
grid on
xlabel('Real'); ylabel('Imag');
