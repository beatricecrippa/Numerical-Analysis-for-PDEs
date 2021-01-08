 [N,Ki,Knn,Kbnn]=textread('conds_schurPN_new','%f %f %f %f','headerlines',2);
 Ki./N
 Knn./(1+log(N)).^2
 Kbnn./(1+log(N)).^2
pause

[H,Ki,Knn,Kbnn]=textread('conds_schurPH_new','%f %f %f %f','headerlines',2);
Ki.*H.^2
Knn.*H.^2
Kbnn

