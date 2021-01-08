function [errors,solutions,femregion,Data]= main2D(TestName,nRef)
%==========================================================================
% Solution of the Poisson's problem with DG finite elements
% (non homogeneous Dirichlet boundary conditions)
%==========================================================================
%
%    INPUT:
%          Data        : (struct)  see dati.m
%          nRef        : (int)     refinement level
%
%    OUTPUT:
%          errors      : (struct) contains the computed errors
%          solutions   : (sparse) nodal values of the computed and exact
%                        solution
%          femregion   : (struct) infos about finite elements
%                        discretization
%          Data        : (struct)  see dati.m
%          
% Usage: 
%    [errors,solutions,femregion,Dati] = main2D('Test1',3)


addpath Assembly
addpath Errors
addpath MeshGeneration
addpath FESpace
addpath PostProcessing


%==========================================================================
% LOAD DATA FOR TEST CASE
%==========================================================================

Data = dati(TestName);

%==========================================================================
% MESH GENERATION
%==========================================================================

[region] = generate_mesh(Data,nRef);

%==========================================================================
% FINITE ELEMENT REGION
%==========================================================================

[femregion] = create_dof(Data,region);


%==========================================================================
% CONNECTIVITY FOR NEIGHBOURING ELEMENTS
%==========================================================================

[neighbour] = neighbours(femregion);


%==========================================================================
% BUILD FINITE ELEMENT MATRICES and RIGHT-HAND SIDE
%==========================================================================

[Matrices] = matrix2D(femregion,neighbour,Data);

%==========================================================================
% SOLVE THE LINEAR SYSTEM
%==========================================================================

u_h = Matrices.A\Matrices.f;


%==========================================================================
% POST-PROCESSING OF THE SOLUTION
%==========================================================================

[solutions]= postprocessing(femregion,Data,u_h);

%==========================================================================
% ERROR ANALYSIS
%==========================================================================
[errors]= compute_errors(Data,femregion,solutions,Matrices.S);


%==========================================================================
% CONDITION NUMBER OF A
%==========================================================================
Data.condA = condest(Matrices.A);







