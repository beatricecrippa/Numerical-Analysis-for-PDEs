function [uex,uexx,uexy,uexz,ff,gam]=setfun_lap_3d
%  SETFUN_LAP_3D  Sets functions and coefficients for lap_3d
%
%    [uex,uexx,uexy,uexz,ff,gam]=setfun_lap_3d
%
% Output: uex =@(x,y,z)[....] function handle to the expression of exact
%               solution
%         uexx =@(x,y,z)[....] function handle to the expression of the first 
%               x-derivative of the exact solution
%         uexy =@(x,y,z)[....] function handle to the expression of the first 
%               y-derivative of the exact solution
%         uexz =@(x,y,z)[....] function handle to the expression of the first 
%               z-derivative of the exact solution
%         ff =@(x,y,z)[....] function handle to the expression of function 
%              at right hand side
%         gam = coefficient of zero-order term (constant >=0)
%

%   Written by Paola Gervasio
%   $Date: 2007/04/01$


syms  x y z
gam=0;
Uex=(x.^2+y.^3+3*z);% sin(2*pi*x)*cos(2*pi*y)*sin(2*pi*z); % exact solution
Uexx=diff(Uex,x); % first x-derivative of exact solution
Uexy=diff(Uex,y); % first y-derivative of exact solution
Uexz=diff(Uex,z); % first z-derivative of exact solution
Uexx2=diff(Uexx,x); % second x-derivative of exact solution
Uexy2=diff(Uexy,y); % second y-derivative of exact solution
Uexz2=diff(Uexz,z); % second z-derivative of exact solution
Ff=(gam*(Uex)-(Uexx2+Uexy2+Uexz2)); % right hand side

uex=strrep(char(Uex),'*','.*');
uexx=strrep(char(Uexx),'*','.*');
uexy=strrep(char(Uexy),'*','.*');
uexz=strrep(char(Uexz),'*','.*');
ff=strrep(char(Ff),'*','.*');

uex=strrep(char(uex),'/','./');
uexx=strrep(char(uexx),'/','./');
uexy=strrep(char(uexy),'/','./');
uexz=strrep(char(uexz),'/','./');
ff=strrep(char(ff),'/','./');

uex=strrep(char(uex),'^','.^');
uexx=strrep(char(uexx),'^','.^');
uexy=strrep(char(uexy),'^','.^');
uexz=strrep(char(uexz),'^','.^');
ff=strrep(char(ff),'^','.^');


uex=@(x,y,z)[eval(uex)];
uexx=@(x,y,z)[eval(uexx)];
uexy=@(x,y,z)[eval(uexy)];
uexz=@(x,y,z)[eval(uexz)];
ff=@(x,y,z)[eval(ff)];

