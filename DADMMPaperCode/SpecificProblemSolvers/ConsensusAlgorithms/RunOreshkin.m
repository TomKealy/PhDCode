% Script that runs Oreshkin for the average consensus:
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
%direct_Oreshkin = '../';
% =========================================================================

% =========================================================================
% Selecting the network

cd(direct_networks);
load Nets_50_nodes.mat;     % File with randomly generated networks
cd(direct_current);

net_num = 1;
Adj = Networks{net_num}.Adj;                   % Adjacency matrix

P = length(Adj);                               % Number of nodes

% Construct the cell neighbors, where neighbors{p} is a vector with the
% neighbors of node p
neighbors = cell(P,1);
Dp = zeros(P,1);

for p = 1 : P
    neighbors{p} = find(Adj(p,:));
    Dp(p) = length(neighbors{p})+1;
end

for p = 1 : P
    neighbors{p} = find(Adj(p,:));
    Dp(p) = length(neighbors{p});
end


% Create struct with network information
vars_network = struct('P', {P}, ...
    'neighbors', {neighbors} ...
    );
% =========================================================================


% =========================================================================
% Selecting the data

% We will use one data point per node
% theta = 100*randn(P,1) + 10*ones(P,1);
% 
% solution = sum(theta)/P;
% 
cd(direct_data);
load Data_Av_Consensus_50_nodes.mat;
cd(direct_current);

% ============================================
% Construct Metropolis-Hastings Weight Matrix
W = zeros(P,P);

for p = 1 : P
    neighbs = neighbors{p};
    Dp_p = Dp(p);
    quant = 0;
    for j_Dp = 1 : Dp_p
        neighbor = neighbs(j_Dp);
        max_deg = max(Dp(p), Dp(neighbor));
        W(p,neighbor) = 1/(max_deg + 1);
        W(neighbor,p) = W(p,neighbor);
        quant = quant + 1/(max_deg + 1);
    end
    W(p,p) = 1 - quant;
end
% ============================================

eigen_vals = eig(W);
l2 = eigen_vals(P-1);

par1 = -0.5;
par2 = 0;
par3 = 3/2;

weight = -(1/(par2 + (par3-1)*l2)^2)*( (par3-1)*l2^2 + par2*l2 + 2*par1 ...
    + 2*sqrt( par1^2 + par1*l2*(par2 + (par3-1)*l2) ) );


% Create struct with problem data
vars_prob = struct('theta', {theta}, ...
    'par1', {par1}, ...
    'par2', {par2}, ...
    'par3', {par3}, ...
    'weight', {weight}, ...
    'W', {W} ...
    );
% =========================================================================

% =========================================================================
% Execute Oreshkin

% Optional input
 ops = struct('max_iter', {300}, ...
     'x_opt', {solution}, ...
     'eps_opt', {1e-5}, ...
     'error_fun', @RelativeError_Av_Consensus ...
 );

[X, vars_prob, ops_out] = Oreshkin(1, vars_prob, vars_network, ops);
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

