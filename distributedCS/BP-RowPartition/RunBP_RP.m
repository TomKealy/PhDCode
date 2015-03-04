% Script that runs D-ADMM for the Basis Pursuit with row partition:
%
%                minimize    ||x||_1                (1)
%                subject to  Ax = b
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
% in the file Nets_50_nodes.mat;

% =========================================================================
% Directories

clear all;
direct_current = pwd;                                      % Current
direct_networks = '/home/tk12098/Documents/MATLAB/DADMM_paper_Code/GenerateData/Networks';           % Networks
% Compressed Sensing Data
direct_data = '/home/tk12098/Documents/MATLAB/DADMM_paper_Code/GenerateData/ProblemData/CompressedSensing';
direct_DADMM = '/home/tk12098/Documents/MATLAB/DADMM_paper_Code/D-ADMM';                                      % D-ADMM
% =========================================================================

max_iter = 500;
rho = 0.1;

% =========================================================================
% Selecting the network

cd(direct_networks);
load Nets_50_nodes.mat;     % File with randomly generated networks
cd(direct_current);

% ================================
% Geometric network
net_num = 4;
% ================================

Adj = Networks{net_num}.Adj;                   % Adjacency matrix
partition_colors = Networks{net_num}.Partition;% Color partition of network
P = length(Adj);                               % Number of nodes
neighbors = Networks{net_num}.Neighbors;

Network = gengeonet(50, 0.25);


% Create struct with network information
vars_network = struct('P', {P}, ...
    'neighbors', {neighbors}, ...
    'partition_colors', {partition_colors} ...
    );



% =========================================================================


% =========================================================================
% Selecting the data

% We will use the Gaussian data provided in GaussianData.mat. We could
% generate our own data and even use the Sparco toolbox:
% http://www.cs.ubc.ca/labs/scl/sparco/

% cd(direct_data);
% %   load Id7_BPDN;
% %   [m,q] = size(A_BP);
% %   opts = spgSetParms('verbosity',0);
% %   solution = spg_bp(A_BP, b, opts);
% %   fprintf('||A_BP*solution-b|| = %E\n', norm(A_BP*solution-b));
% %   fprintf('norm(solution,1) = %E\n', norm(solution,1));
% load GaussianData.mat;             % This contains:
%                                    % A_BP      - the matrix A in (1)
%                                    % b         - vector b in (1)
%                                    % k         - the sparsity level
%                                    % m         - number of rows of A
%                                    % q         - number of columns of A
%                                    % solution  - solution of (1)
% cd(direct_current);

L=201; %number of sub-bands
L0 = (L-1)/2;
m=50;
positions = randi(L,[1,5]);%generate random spikes for signal

Tx_psd = zeros(1,L); %Tx PSD
Tx_psd(positions) = 1;
% figure
% plot(Tx_psd)

Eb_N0_dB = [-10:1:10];
[r,c] = size(Eb_N0_dB);
mse_admm_av = zeros(1,c);
mse_soln_av = zeros(1,c);
mse_gpsr_av = zeros(1,c);
S = randn(50,L);
%for j=1:c

runs = 10
mse_admm = zeros(1,runs);
mse_soln = zeros(1,runs);
mse_gpsr = zeros(1,runs);

for i=1:runs
    
    h = (1/sqrt(2))*(randn(1,L)+1i*randn(1,L));
    
    H = diag(h);
    
    n = 1/sqrt(2)*[randn(1,m) + j*randn(1,m)]; % white gaussian noise, 0dB variance
    
    theta = exp(-1i*2*pi/L);
    F = theta.^([0:L-1]'*[-L0:L0]);
    np = 1:L0;
    nn = (-L0):1:-1;
    % This is for digital input only. Note that when R -> infinity,
    % D then coincides with that of the paper
    dn = [   (1-theta.^nn)./(1-theta.^(nn))/(L)      1/L    (1-theta.^np)./(1-theta.^(np))/(L)];
    D = diag(dn);
    
    A = S;
    b = A*Tx_psd';
    
    central_solution = spg_bp(A, b);
    beta = 0.3
    solution = GPSR_BB(b, A, beta, 'Verbose', 0, 'StopCriterion', 3, ...
        'ToleranceA', 1e-2);
    
    mse_soln(i) = norm(central_solution' - Tx_psd)/norm(Tx_psd)
    mse_gpsr(i) = norm(solution' - Tx_psd)/norm(Tx_psd)
    %Check if matrix partition is possible (all blocks with same size)
    
    %ADMM Soln
    if mod(m,P) ~= 0
        error('m divided by P must be integer');
    end
    m_p = P/P;                         % Number of rows of A each node stores
    
    %The function (that should be built by the user of D-ADMM) %
    %minimize_quad_prog_plus_l1_BB solves the problem
    %min {(1/P)||x||_1 + v'*x + c||x||^2 : Ap x = bp}.
    %We use duality to solve that problem and we need to initialize the dual
    %variable before using the function. There is one dual variable per node.
    
    dual_var = cell(P,1);
    
    for p = 1 : P                      % Initialization with zeros
        dual_var{p} = zeros(m_p,1);
    end
    
    %Create struct with problem data used in 'minimize_quad_prog_plus_l1_BB'
    vars_prob = struct('handler', @minimize_quad_prog_plus_l1_BB,...
        'A_BP', {A}, ...
        'b_BP', {DigitalSamples},...
        'dual_var', {dual_var} ...
        );
    
    ops = struct('rho', {rho}, ...
        'max_iter', {max_iter}, ...
        'eps', {1e-4}, ...
        'x_opt', {u},...
        'eps_opt', {1e-4}, ...
        'turn_off_eps', {0} ...
        );
    
    cd(direct_DADMM);
    [X, vars_prob, ops_out] = DADMM(L, vars_prob, vars_network, ops);
    cd(direct_current);
    
    iterations = ops_out.iterations;
    stop_crit = ops_out.stop_crit;
    error_iterations = ops_out.error_iterations;
    iter_for_errors = ops_out.iter_for_errors;
    
    fprintf('norm(A*X{1} - b)/norm(b) = %f\n',norm(A*X{1} - b)/norm(b));
    fprintf('||X{1}||_1 = %f\n',norm(X{1},1));
    fprintf('||X{1} - solutionl||/||solution|| = %f\n',norm(X{1}-solution)/norm(solution));
    fprintf('Number of iterations = %d\n', iterations);
    fprintf('stop_crit = %s\n', stop_crit);
    fprintf('iter_for_errors = \n');
    num_rows = size(iter_for_errors, 1);
    for i_g = 1 : num_rows
        fprintf('%E    %d\n', iter_for_errors(i_g,1), iter_for_errors(i_g,2));
    end
    
    mse_admm(i) = norm(X{1}' - Tx_psd)/norm(Tx_psd)
    
    
    figure(1);clf;
    semilogy(1:iterations,error_iterations(1:iterations), 'b');
    title('error\_{iterations}');
    
    figure(2)
    plot(solution)
    
    figure(3)
    plot(X{1})
end
%     =========================================================================
mse_admm_av(j) = sum(mse_admm)/runs;
mse_soln_av(j) = sum(mse_soln)/runs;
mse_gpsr_av(j) = sum(mse_soln)/runs;
%end

figure
plot(Eb_N0_dB, mse_soln_av, '-ro', Eb_N0_dB, mse_admm_av,'-xb')
hleg = legend('cetralised','distributed')
xlabel('Eb-N0-dB')
ylabel('mse')

figure
plot(X{1})
xlabel('f')
ylabel('PSD')

figure
plot(central_solution)
xlabel('f')
ylabel('PSD')

