# Numerical-Analysis-for-PDEs

Laplace problems solved via neural networks.

In problems.py four different classes of Laplace problems are defined:
- 2D with 0 Dirichlet boundary data
- 2D without boundary data
- 2D with varying Dirichlet boundary data
- N-dimensional with varying boundary data

In pdebase.py classes for the solution of the previous problems via neural networks are defined:
- 2D neural network with 0 Dirichlet boundary data
- 2D neural network with varying Dirichlet boundary data
- N-dimensional neural network with varying boundary data.
One subnetwork works on the boundary data and one on the inner domain: the subnetwokrs are dense neural networks with tanh activation and output dimensionality 256; at each neuron the sum of square error (SSE) is minimised.

Run problem1.py, problem2.py, problem3.py and problemN.py to get the solution and display the approxiamtion errors.
Before running the problems, create in your directory folders named "p1", "p2", "p3", "pN" respleively, where the graphs will be saved.
