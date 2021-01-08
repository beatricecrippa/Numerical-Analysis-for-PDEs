% FIG4_8  Script to produce Fig 4.8, pag. 201 CHQZ2
%
% Chebyshev collocation first-derivative eigenvalues
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
    npdx=nx+1;
    [x,wx]=xwcgl(npdx); dx=dercgl(x,npdx);
A=dx(1:nx,1:nx);
d=eig(A);
fig=figure(...,
    'Name',['Fig. 4.8 N=',num2str(nx)],...
    'Visible','on');
plot(real(d),imag(d),'ko'); grid on
xlabel('Real'); ylabel('Imag');
axis equal
legend(['N=',num2str(nx)]);
end
