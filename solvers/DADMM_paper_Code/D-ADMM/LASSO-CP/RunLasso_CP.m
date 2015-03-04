% Script that runs D-ADMM for LASSO with row partition:
%
%          minimize    ||x||_1                  (1)
%          subject to  ||Ax - b|| <= sigma
%
% where x is the variable and ||.||_1 is the L1 norm. Sigma and b are known
% by all nodes, and the matrix A in partitioned by columns:
%
%          A = [A1 A2 ... AP]
%
% and distributed across the P nodes in the network. 

% =========================================================================
% Directories

direct_current = pwd;                                      % Current
direct_networks = '../../GenerateData/Networks';           % Networks
% Compressed Sensing Data
direct_data = '../../GenerateData/ProblemData/CompressedSensing';
direct_DADMM = '../';                                      % D-ADMM
% =========================================================================


% =========================================================================
% Selecting the network

cd(direct_networks);
load Nets_50_nodes.mat;     % File with randomly generated networks
cd(direct_current);

% Select the network number: 1 to 7
net_num = 1;

Adj = Networks{net_num}.Adj;                   % Adjacency matrix
partition_colors = Networks{net_num}.Partition;% Color partition of network
neighbors = Networks{net_num}.Neighbors;

P = length(Adj);                               % Number of nodes

% Create struct with network information
vars_network = struct('P', {P}, ...
    'neighbors', {neighbors}, ...
    'partition_colors', {partition_colors} ...
    );
% =========================================================================


% =========================================================================
% Selecting the data

% We use the Sparco toolbox: http://www.cs.ubc.ca/labs/scl/sparco/ or our
% own generated data (Id = 0)

Id = 902;

delta = 1e-2;                                % Regularization factor
sigma = 1e-1;                                 % Parameter of LASSO

cd(direct_data);
if Id == 0        % Gaussian data
    load GaussianData.mat;
    [m, n] = size(A_BP);
else
    Prob = generateProblem(Id);
    A_BP_aux = classOp(Prob.A);         % (Almost) explicit matrix A_BP
    b = Prob.b;                         % Vector b
    [m, n] = size(A_BP_aux);
    
    % Get the explicit matrix A_BP (in double format)
    A_BP = zeros(m, n);
    for i = 1 : m
        ei = zeros(m,1);
        ei(i) = 1;
        A_BP(i,:) = (A_BP_aux'*ei)';
    end
    clear A_BP_aux;
end
cd(direct_current);

% Scale b (so that x has large components) and noise to b
b = 100*b + 0.1*sigma*randn(m,1);
                                   
% For the groundtruth, we use the spgl1 solver
options_spgl1 = struct('verbosity', {0});
solution = spgl1(A_BP, b, 0, sigma, [], options_spgl1);

fprintf('||A_BP*solution-b|| = %E < %E\n', norm(A_BP*solution-b), sigma);
fprintf('norm(solution,1) = %E\n', norm(solution,1));
fprintf('card(solution) = %d\n', sum(solution~=0));


% Check if matrix partition is possible (all blocks with same size)
if mod(n,P) ~= 0
    error('q divided by P must be integer');
end
n_p = n/P;                         % Number of rows of A each node stores

A = cell(P,1);

% =========================================================================
% Create structure with data

max_sing_val_sq = cell(P,1);
previous_t = cell(P,1);
x_opt = cell(P,1);                 % Estimates of x along the algorithm

for p = 1 : P
    A{p} = A_BP(: , (p-1)*n_p + 1 : p*n_p );
    sing_vals = svd(A{p});
    max_sing_val_sq{p} = max(sing_vals)^2;
    previous_t{p} = 0;
    x_opt{p} = zeros(n_p,1);
end

one_over_P_b = 1/P*b;
sig_P = sigma/P;

% Create struct with problem data used in 'minimize_quad_prog_plus_l1_BB'
vars_prob = struct('handler', @NesterovForLassoSubProb,...
    'n', {n_p}, ...
    'A', {A}, ...
    'max_sing_val_sq', {max_sing_val_sq}, ...
    'delta', {delta}, ...
    'one_over_P_b', {one_over_P_b}, ...
    'sig_P', {sig_P}, ...
    'x_opt', {x_opt}, ...
    'previous_t', {previous_t} ...
    );
% =========================================================================


% =========================================================================
% Execute D-ADMM

% Optional input
 ops = struct('rho', {1}, ...
     'max_iter', {1000}, ...
     'x_opt', {solution}, ...
     'eps_opt', {2e-3}, ...
     'error_fun', @error_function_Lasso_for_DADMM ...
 );

cd(direct_DADMM);
[X, vars_prob, ops_out] = DADMM(m, vars_prob, vars_network, ops);
cd(direct_current);
% =========================================================================

% =========================================================================
% Print results

iterations = ops_out.iterations;
stop_crit = ops_out.stop_crit;
error_iterations = ops_out.error_iterations;
iter_for_errors = ops_out.iter_for_errors;

% Recover solution
n = P*n_p;
est = vars_prob.x_opt;
x_estimate = zeros(n,1);
for p = 1 : P 
    x_estimate( (p-1)*n_p + 1 : p*n_p ) = est{p};
end

fprintf('norm(A*x_estimate - b)/norm(b) = %f\n',norm(A_BP*x_estimate - b)/norm(b));
fprintf('||x_estimate||_1 = %f\n',norm(x_estimate,1));
fprintf('||x_estimate - solutionl||/||solution|| = %f\n',norm(x_estimate-solution)/norm(solution));
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






