function [A,f]=C_matrix2D(Dati,femregion)
%% [A,f] = C_matrix2D(Dati,femregion)
%==========================================================================
% Assembly of the stiffness matrix A and rhs f
%==========================================================================
%    called in C_main2D.m
%
%    INPUT:
%          Dati        : (struct)  see C_dati.m
%          femregion   : (struct)  see C_create_femregion.m
%
%    OUTPUT:
%          A           : (sparse(ndof,ndof) real) stiffnes matrix
%          f           : (sparse(ndof,1) real) rhs vector


addpath FESpace
addpath Assembly

fprintf('============================================================\n')
fprintf('Assembling matrices and right hand side ... \n');
fprintf('============================================================\n')


% connectivity infos
ndof         = femregion.ndof; % degrees of freedom
nln          = femregion.nln;  % local degrees of freedom
ne           = femregion.ne;   % number of elements
connectivity = femregion.connectivity; % connectivity matrix


% shape functions
[basis] = C_shape_basis(Dati);

% quadrature nodes and weights for integrals
[nodes_2D, w_2D] = C_quadrature(Dati);

% evaluation of shape bases 
[dphiq,Grad] = C_evalshape(basis,nodes_2D);


% Assembly begin ...
A = sparse(ndof,ndof);  % Global Stiffness matrix
f = sparse(ndof,1);     % Global Load vector

for ie = 1 : ne
     
    % Local to global map --> To be used in the assembly phase
    iglo = connectivity(1:nln,ie);
  
    
    [BJ, pphys_2D] = C_get_Jacobian(femregion.coord(iglo,:), nodes_2D, Dati.MeshType);
    % BJ        = Jacobian of the elemental map 
    % pphys_2D = vertex coordinates in the physical domain 
   
    %=============================================================%
    % STIFFNESS MATRIX
    %=============================================================%
    
    % Local stiffness matrix 
    [A_loc] = C_lap_loc(Grad,w_2D,nln,BJ);
    [M_loc] = C_mass_loc(dphiq,w_2D,nln,BJ);
    [Adv_loc]=C_adv_loc(Grad,dphiq, Dati.beta,w_2D,nln,BJ);

    % Assembly phase for stiffness matrix
    A(iglo,iglo) = A(iglo,iglo) + Dati.mu*A_loc + Dati.sigma*M_loc +Adv_loc;
    
    %==============================================
    % FORCING TERM --RHS
    %==============================================

    % Local load vector
    [load] = C_loc_rhs2D(Dati.force,dphiq,BJ,w_2D,pphys_2D,nln);    

    % Assembly phase for the load vector
    f(iglo) = f(iglo) + load;

end
