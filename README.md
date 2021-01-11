# Numerical-Analysis-for-PDEs

1. Python folder:
Laplace problems solved via neural networks.

In pdebase.py base classes of neural networks with methods for the solution of the Laplace problem are defined:
- 2D neural network with Dirichlet boundary condition approximated by a neural network on the boundary (NNPDE2)
- 2D neural network with fixed Dirichlet boundary condition (NNPDE)
- N-dimensional neural network with Dirichlet boundary condition approximated by a neural network on the boundary (NNPDE_ND).

In problems.py classes of specific problems, derived from pdebase.py classes, are defined:
- Smooth solution
- Solution with a peak
- Less regular solution

in all the three previous cases.

In particular,
- Problem_1, ProblemPeak and ProblemBLSingularity describe the previous problems as derived from NNPDE
- Problem_1_BD, ProblemPeak_BD and ProblemBLSingularity_BD describe the previous problems as derived from NNPDE2
- HighDimensionSmooth, HighDimensionPeak, HighDimensionSingularity describe the previous problems as derived from NNPDE_ND (N-dimensional version).

One subnetwork works on the boundary data and one on the inner domain: the subnetwokrs are dense neural networks with tanh activation and output dimensionality 256; at each neuron the sum of square error (SSE) is minimised.

Run problem1.py, problem2.py, problem3.py and problemN1.py, problemN2.py, problemN3.py to get the solutions and display the approxiamtion errors.
Before running the problems, create in your directory folders named "p1", "p2", "p3", "high/p1", "high/p2" and "high/p3" respleively, where the graphs will be saved.

2. Matlab folder:
Laplace problems solved via traditional numerical method.

In the folder CG_FEM the Galerkin Finite Elements Method is applied to the three previous problems.
Run C_Convergence_test with parameters 'Test1', 'Test2', 'Test3' respectively for the problem with smooth solution, solution with peak and solution with singularity.

