This folder contains the file DLassoCompos.m which implements the version of D-Lasso that appears in

[1]  G. Mateos, J. Bazerque, and G. Giannakis, "Distributed Sparse Linear Regression," IEEE Trans. Sig. Proc., Vol. 58, No. 10, 2010

The function DLassoCompos.m implements Algorithm 3 of [1]. Although [1] applies the algorithm to a specific problem, the algorithm can be applied to solve 

    minimize \sum_{p=1}^P  g_p(x) + h_p(x)
       x

in a distributed way, in a network with P nodes. The functions g_p and h_p are only known to node p. The rationale behind the method is to clone the variable x throughout the nodes (copy x_p is assigned to node p) and, at the same time, clone the variable x_p into y_p, i.e., the cloning of x is associated to the network, and the cloning of x_p into y_p is associated to node p only. After this cloning, we derive the algorithm as in the regular D-Lasso, but we also dualize the constraints x_p = y_p: the dual variable associated to x_p = y_p is updated only by node p. 

[1] applies this problem to h_p(x) = 0.5*||Ap*x - bp||^2 (will be associated with y),  g_p(x) = 0.5*beta/P||x||_1 (will be associated with x).

