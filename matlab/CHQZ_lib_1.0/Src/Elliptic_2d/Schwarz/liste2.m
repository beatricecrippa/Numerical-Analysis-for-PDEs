function [listaint,listadir]=liste2(ifro);
% LISTE2 Similar to ~/Lelve_2/liste.m 
% 
% called by stiffq1.m

%   Written by Paola Gervasio
%   $Date: 2007/04/01$


noe=length(ifro);
nint=0;
listaint=[];listadir=[];
for i=1:noe
if(ifro(i)~=1)
nint=nint+1;
listaint=[listaint;i];
else
listadir=[listadir;i];
end
end
