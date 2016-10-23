% Script that runs D-ADMM for the average consensus:
%
%          minimize    sum_{p=1}^P (x - theta_p)^2               (1)
%             x
%
% We use a network provided in the file Nets_50_nodes.mat;

% =========================================================================
% Directories

direct_current = pwd;                                      % Current
direct_networks = '../../GenerateData/Networks';           % Networks
direct_data = '../../GenerateData/ProblemData/Consensus';  % Consensus Data
direct_DADMM = '../';                                      % D-ADMM
% =========================================================================


% =========================================================================
% Selecting the network

cd(direct_networks);
load Nets_50_nodes.mat;             % File with randomly generated networks
cd(direct_current);

net_num = 1;

Adj = Networks{net_num}.Adj;                   % Adjacency matrix
partition_colors = Networks{net_num}.Partition;% Color partition of network
P = length(Adj);                               % Number of nodes
neighbors = Networks{net_num}.Neighbors;
optim_rho = Networks{net_num}.OptimRho_Cons;
%optim_rho = 1;

% Create struct with network information
vars_network = struct('P', {P}, ...
    'neighbors', {neighbors}, ...
    'partition_colors', {partition_colors} ...
    );
% =========================================================================


% =========================================================================
% Selecting the data

% We will use one data point per node
% theta = 100*randn(P,1) + 10*ones(P,1);
% 
% solution = sum(theta)/P;

cd(direct_data);
load Data_Av_Consensus_50_nodes.mat;
cd(direct_current);
                                   
% Create struct with problem data
vars_prob = struct('handler', @AverageConsensusSubProb,...
    'theta', {theta} ...
    );
% =========================================================================

% =========================================================================
% Execute D-ADMM

% Optional input
 ops = struct('rho', {optim_rho}, ...
     'max_iter', {300}, ...
     'x_opt', {solution}, ...
     'eps_opt', {1e-5}, ...
     'error_fun', @RelativeError_Av_Consensus ...
 );

cd(direct_DADMM);
[X, vars_prob, ops_out] = DADMM(1, vars_prob, vars_network, ops);
cd(direct_current);
% =========================================================================

% =========================================================================
% Print results

iterations = ops_out.iterations;
stop_crit = ops_out.stop_crit;
error_iterations = ops_out.error_iterations;
iter_for_errors = ops_out.iter_for_errors;

% Recover solution
x_estimate = zeros(P,1);
for p = 1 : P
    x_estimate(p) = X{p};
end

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






