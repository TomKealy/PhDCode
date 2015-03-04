% Script that runs D-Lasso-Compos for BPDN with row partition:
%
%          minimize    0.5*||A*x - b||^2 + beta*||x||_1           (1)
%             x
%
% where x is the variable and ||.||_1 is the L1 norm. Row partition
% (VP) means that the matrix A and the vector b are partitioned as:
%
%               [    A1    ]              [b1]
%               [    A2    ]              [b2]
%          A =  [   ...    ]          b = [..]
%               [    AP    ]              [bP]
%
% and node p in the network only knows Ap and bp. We use a network provided
% in the file Nets_10_nodes.mat;


% =========================================================================
% Directories

direct_current = pwd;                                      % Current
direct_networks = '../../../GenerateData/Networks';           % Networks
% Compressed Sensing Data
direct_data = '../../../GenerateData/ProblemData/CompressedSensing';
direct_DLassoComp = '../';                                 % D-LassoComp
% =========================================================================

rho = 0.1;
MAX_ITER = 1000;
EPS_OPT = 1e-4;

% =========================================================================
% Selecting the network

cd(direct_networks);
load Nets_10_nodes.mat;     % File with randomly generated networks
cd(direct_current);

% Select the network number: 1 to 7
net_num = 1;

Adj = Networks{net_num}.Adj;                   % Adjacency matrix
P = length(Adj);                               % Number of nodes

% Construct the cell neighbors, where neighbors{p} is a vector with the
% neighbors of node p
neighbors = cell(P,1);

for p = 1 : P
    neighbors{p} = find(Adj(p,:));
end

% Create struct with network information
vars_network = struct('P', {P}, ...
    'neighbors', {neighbors} ...    
    );
% =========================================================================


% =========================================================================
% Selecting the data

% We use the Sparco toolbox: http://www.cs.ubc.ca/labs/scl/sparco/ or our
% own generated data (Id = 0)

Id = 902;

beta = 0.3;                                   

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

% Add Gaussian, zero-mean noise to vector b
std_noise = 0.5*norm(b);
b = b + std_noise*randn(m,1);


% For the groundtruth, we use the GPSR solver
[solution] = GPSR_BB(b, A_BP, beta, 'Verbose', 0, 'StopCriterion', 3, ...
    'ToleranceA', 1e-6);

fprintf('||A_BP*solution-b|| = %E\n', norm(A_BP*solution-b));
fprintf('norm(solution,1) = %E\n', norm(solution,1));


% Check if matrix partition is possible (all blocks with same size)
if mod(m,P) ~= 0
    error('m divided by P must be integer');
end
m_p = m/P;                         % Number of rows of A each node stores

bpAp = cell(P,1);       % Will contain Ap'*bp
M_inv = cell(P,1);      % Will contain inv(Ap'*Ap + rho*eye(n))

for p = 1 : P
    Ap = A_BP(1+(p-1)*m_p : p*m_p , :);
    bp = b(1+(p-1)*m_p : p*m_p);
    bpAp{p} = Ap'*bp;
    M_inv{p} = inv(Ap'*Ap + rho*eye(n));
end
                                   
vars_prob = struct('handler_cons', @Solver_BPDN_g, ...
    'handler_uncons', @Solver_BPDN_h, ...
    'M_inv', {M_inv}, ...
    'Apbp', {bpAp}, ...
    'P', {P}, ...
    'beta', {beta} ...
    );

% =========================================================================

% =========================================================================
% Execute D-Lasso-Compos

% Optional input
 ops = struct('rho', {rho}, ...
     'max_iter', {MAX_ITER}, ...
     'x_opt', {solution}, ...
     'eps_opt', {EPS_OPT} ...
 );

cd(direct_DLassoComp);
fprintf('DLassoCompos: Start\n');
[X, vars_prob, ops_out] = DLassoCompos(n, vars_prob, vars_network, ops);
fprintf('DLassoCompos: Finish\n');
cd(direct_current);
% =========================================================================

% =========================================================================
% Print results

iterations = ops_out.iterations;
stop_crit = ops_out.stop_crit;
error_iterations = ops_out.error_iterations;
iter_for_errors = ops_out.iter_for_errors;

fprintf('norm(A*X{1} - b)/norm(b) = %f\n',norm(A_BP*X{1} - b)/norm(b));
fprintf('||X{1}||_1 = %f\n',norm(X{1},1));
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






