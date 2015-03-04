This folder contains the following files:

   RunBP_RP.m - A script that uses D-ADMM to solve Basis Pursuit with 
                row partition. It uses networks from the file 
                Nets_50_nodes.mat. By default, it does not require
                external Matlab software, but it can be set to use data
                from the Sparco toolbox:
                http://www.cs.ubc.ca/labs/scl/sparco/                   
   
   minimize_quad_prog_plus_l1_BB.m - Solver for D-ADMM.
   

The file TEST_ADMM_VS_BB.m compares the Barzilai-Borwein algorithm with ADMM for computing the
prox of ||x||_1 + linear, subject to linear constraints. BB wins although not by much.

