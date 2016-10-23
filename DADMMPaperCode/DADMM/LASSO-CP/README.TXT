This folder contains the following files:

   RunLasso_CP.m - A script that uses D-ADMM to solve LASSO with column
                   partition. It uses networks from the file 
                   Nets_50_nodes.mat. By default, it requires the Sparco 
                   toolbox to be installed. See: 
                   http://www.cs.ubc.ca/labs/scl/sparco/ 
                   For computing the optimal solution, this script also uses the
                   spgl1 solver. See:
                   http://www.cs.ubc.ca/~mpf/spgl1/                  
   
   NesterovForLassoSubProb.m - Solver for D-ADMM. It uses a Nesterov algorithm.

   error_function_Lasso_for_DADMM.m - Computes the relative error at each 
                                      D-ADMM iteration.


