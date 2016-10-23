% Script that runs D-ADMM for the following Support Vector Machine (SVM)
% formulation: given a set of K pairs of points: (x_k,y_k), where x_k is a 
% vector n x 1 and y_k is a scalar {-1,+1}, find a linear separator 
% {(s,r): y = s'*x + r} that ``best'' separates the points x_k associated 
% to y_k = 1 and y_k = -1. By best, we mean the solution of
%
%                minimize    0.5*||s||^2 + beta*ones(1,K)*z        (1)
%                subject to  y_k*(s'*x_k - r) >= 1 - z_k,  k = 1,...,K
%                            z >= 0
%
% where the variable is (s, r, z), and beta is a constant. To apply it to 
% D-ADMM, we consider the data is distributed among P entities. Namely, 
% some of the points (x_k,y_k) are with the pth entity.

% =========================================================================
% Directories

direct_current = pwd;                                      % Current
direct_networks = '../../GenerateData/Networks';           % Networks
direct_data = '../../GenerateData/ProblemData/SVMs';       % SVMs Data
direct_DADMM = '../';                                      % D-ADMM
% =========================================================================


% =========================================================================
% Selecting the network

cd(direct_networks);
load Nets_50_nodes.mat;     % File with randomly generated networks
cd(direct_current);

rho = 0.5;
MAX_ITER = 5000;
eps_opt = 1e-3;

% Select the network number: 1 to 7
net_num = 1;

Adj = Networks{net_num}.Adj;                   % Adjacency matrix
partition_colors = Networks{net_num}.Partition;% Color partition of network

P = length(Adj);                               % Number of nodes

% Construct the cell neighbors, where neighbors{p} is a vector with the
% neighbors of node p
neighbors = cell(P,1);

Dp = zeros(P,1);
for p = 1 : P
    neighbors{p} = find(Adj(p,:));
    Dp(p) = length(neighbors{p});
end

% Create struct with network information
vars_network = struct('P', {P}, ...
    'neighbors', {neighbors}, ...
    'partition_colors', {partition_colors} ...
    );
% =========================================================================


% =========================================================================
% Selecting the data

% We use the database from the UCI Machine Learning Repository
% http://archive.ics.uci.edu/ml/

beta = 1;

cd(direct_data);
load Iris_overlap.mat;
cd(direct_current);

X = A_mat;
Y_vec = D_vec;
[K , n] = size(A_mat);


% Check if all nodes can store the same ammount of data
if mod(K,P) ~= 0
    error('K divided by P must be integer');
end
m = K/P;                         % Number of rows of A each node stores

% =================================================
% Compute the solution with quadprog
YX = diag(Y_vec)*X;
Y1 = Y_vec;
Id_K = eye(K);

A_qp = [-YX , Y1 , -Id_K ; zeros(K,n+1) , -Id_K];
b_qp = [-ones(K,1) ; zeros(K,1)];
H_qp = zeros(n+1+K);
H_qp(1:n,1:n) = eye(n);
f_qp = [zeros(n+1,1) ; beta*ones(K,1)];

x_qp = quadprog(H_qp,f_qp,A_qp,b_qp);

solution = x_qp(1:n+1);
% =================================================

% =================================================
% Compute the solution with cvx
% cvx_begin
%     variables s_cvx(n) r_cvx xi_cvx(K);
%     minimize(0.5*square_pos(norm(s_cvx)) + beta*ones(1,K)*xi_cvx);
%     subject to
%         YX*s_cvx - Y1*r_cvx >= ones(K,1) - xi_cvx;
%         xi_cvx >= zeros(K,1);
% cvx_end
% 
% solution = [s_cvx ; r_cvx];
% 
% cvx_quiet(true);
% =================================================


% *************************************************************************
% COMMENT THE NEXT CODE, IF YOU WANT THE BARZILAI-BORWEIN SOLVER
% *************************************************************************
A_qp = cell(P,1);
b_qp = [-ones(m,1) ; zeros(m,1)];  % b_p is the same for all nodes
Xi_p = cell(P,1);                  % Holds the slack variables of each node

Id_m = eye(m);
for p = 1 : P
    ind = (1 + (p-1)*m) : p*m;
    YX_p = YX(ind,:);
    Y1_p = Y1(ind,:);
        
    A_qp{p} = [-YX_p , Y1_p , -Id_m ; zeros(m,n+1) , -Id_m];
    Xi_p{p} = zeros(m,1);
end

% Create struct with problem data
vars_prob = struct('handler', @SVM_solver, ...
    'A_qp', {A_qp}, ...                             
    'b_qp', {b_qp}, ...                             
    'Xi_p', {Xi_p}, ...                             
    'beta', {beta}, ...                             
    'm', {m}, ...
    'n', {n}, ...
    'P', {P} ...    
    );
% *************************************************************************
% *************************************************************************

% *************************************************************************
% UNCOMMENT THE NEXT CODE, IF YOU WANT THE BARZILAI-BORWEIN SOLVER
% *************************************************************************
% A_p = cell(P,1);
% warm_z = cell(P,1);
% 
% for p = 1 : P
%     ind = (1 + (p-1)*m) : p*m;
%     Ap = [-YX(ind,:) , Y_vec(ind)];
%     
%     if rank(Ap) == m
%         A_p{p} = Ap;
%     else
%         error('Matrix Ap = [-YX(ind,:) ; Y_vec(ind)] is not full rank');
%     end
%     
%     warm_z{p} = zeros(2*m,1);
%     
% end
% 
%                                                                              
% % Create struct with problem data
% vars_prob = struct('handler', @SVM_solver, ...
%     'beta', {beta}, ...                             
%     'm', {m}, ...
%     'n', {n}, ...
%     'P', {P}, ...
%     'A_p', {A_p}, ...
%     'warm_z', {warm_z} ...
%     );
% *************************************************************************
% *************************************************************************


% =========================================================================

% =========================================================================
% Execute D-ADMM

% Optional input
 ops = struct('rho', {rho}, ...
     'max_iter', {MAX_ITER}, ...
     'x_opt', {solution}, ...
     'eps_opt', {eps_opt} ...     
 );

cd(direct_DADMM);
[X, vars_prob, ops_out] = DADMM(n+1, vars_prob, vars_network, ops);
cd(direct_current);
% =========================================================================

% =========================================================================
% Print results

iterations = ops_out.iterations;
stop_crit = ops_out.stop_crit;
error_iterations = ops_out.error_iterations;
iter_for_errors = ops_out.iter_for_errors;

fprintf('||X{1} - solutionl||/||solution|| = %f\n',norm(X{1}-solution)/norm(solution));
fprintf('Number of iterations = %d\n', iterations);
fprintf('stop_crit = %s\n', stop_crit);
fprintf('iter_for_errors = \n');
num_rows = size(iter_for_errors, 1);
for i_g = 1 : num_rows
    fprintf('%E    %d\n', iter_for_errors(i_g,1), iter_for_errors(i_g,2));
end


figure(1);clf;
semilogy(1:iterations,error_iterations(1:iterations), 'b');
title('error\_{iterations}');
% =========================================================================






