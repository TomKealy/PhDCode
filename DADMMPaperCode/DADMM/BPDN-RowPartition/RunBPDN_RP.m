clear all;
close all;
% Script that runs D-ADMM for BPDN with row partition:
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
% in the file Nets_50_nodes.mat;

% =========================================================================
% Directories

% direct_current = pwd;                                      % Current
% direct_networks = '../../GenerateData/Networks';           % Networks
% % Compressed Sensing Data
% direct_data = '../../GenerateData/ProblemData/CompressedSensing';
% direct_DADMM = '../';                                      % D-ADMM
clear all;
direct_current = pwd;                                      % Current
direct_networks = '/home/tk12098/Documents/MATLAB/DADMM_paper_Code/GenerateData/Networks';           % Networks
% Compressed Sensing Data
direct_data = '/home/tk12098/Documents/MATLAB/DADMM_paper_Code/GenerateData/ProblemData/CompressedSensing';
direct_DADMM = '/home/tk12098/Documents/MATLAB/DADMM_paper_Code/D-ADMM';
% =========================================================================


% =========================================================================
% Selecting the network

cd(direct_networks);
load Nets_50_nodes.mat;     % File with randomly generated networks
cd(direct_current);

% Select the network number: 1 to 7
net_num = 4;

Adj = Networks{net_num}.Adj;                   % Adjacency matrix
partition_colors = Networks{net_num}.Partition;% Color partition of network

P = length(Adj);                               % Number of nodes

% Construct the cell neighbors, where neighbors{p} is a vector with the
% neighbors of node p
neighbors = cell(P,1);

for p = 1 : P
    neighbors{p} = find(Adj(p,:));
end

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

Id = 0;

beta = [0.1, 0.5, 1, 10, 15];
size_beta = size(beta, 2);
    
n=200;
m=50;
L=n;

positions = randi(L,[1,5]);%generate random spikes for signal

Tx_psd = zeros(1,L); %Tx PSD
Tx_psd(positions) = 1;
Eb_N0_dB = 1;
runs = 1;
[r,c] = size(Eb_N0_dB);

mse_admm_av = zeros(1,c);
mse_soln_av = zeros(1,c);

mse_admm = zeros(1,runs);
mse_soln = zeros(1,runs);

max_iter = 1000;

for i=1:size(Eb_N0_dB,1)
    for j = 1:runs
        S = randn(m,L);
        h = sqrt(1/2)*(randn(1,L)+1j*randn(1,L)); %TOK
        H = diag(abs(h)); %TOK
        
        A_BP = S;%*H;
        sigma = 10^(-Eb_N0_dB(i)/20);
        eta = randn(1,m)/m;
        noise_sum = sum(eta);
        b = A_BP*Tx_psd' + sigma*eta';
        
        
        step_size = 0.05 %nthroot(1/(max(abs(eig(A_BP'*A_BP)))),3);
        regression_parameter = 0.1;%
        relaxation_parameter = 1.0;
        % For the groundtruth, we use the spgl1 solver
        
        solution = spgl1(A_BP, b, 0, 0.0001, []);
        
        fprintf('||A_BP*solution-b|| = %E\n', norm(A_BP*solution-b));
        fprintf('norm(solution,1) = %E\n', norm(solution,1));
        
        % Check if matrix partition is possible (all blocks with same size)
        if mod(m,P) ~= 0
            error('m divided by P must be integer');
        end
        m_p = m/P;                         % Number of rows of A each node stores
               
        % Create struct with problem data used in 'minimize_quad_prog_plus_l1_BB'
        vars_prob = struct('A_BPDN', {A_BP}, ...
            'b_BPDN', {b}, ...
            'm_p', {m_p}, ...
            'P', {P}, ...
            'beta', {regression_parameter}, ...
            'relax', {relaxation_parameter} ...
            );
        % =========================================================================
        
        % =========================================================================
        % Execute D-ADMM
        
        % Optional input
        ops = struct('rho', {step_size}, ...
            'max_iter', {max_iter}, ...
            'x_opt', {Tx_psd'}, ...
            'eps_opt', {1e-2}, ...
            'turn_off_eps', {0}....
            );
        
        cd(direct_DADMM);
        tic
        [X, Z, vars_prob, ops_out_accel] = DADMM_accel(n, vars_prob, vars_network, ops);
        toc
        cd(direct_current);
        
        % =========================================================================
        
        % =========================================================================
        % Print results
        
        lambda = sqrt(2*log10(n));
        
         % Create struct with problem data used in 'minimize_quad_prog_plus_l1_BB'
        vars_prob = struct('A_BPDN', {A_BP}, ...
            'b_BPDN', {b}, ...
            'm_p', {m_p}, ...
            'P', {P}, ...
            'beta', {regression_parameter}, ...
            'relax', {relaxation_parameter} ...
            );
        % =========================================================================
        
        % =========================================================================
        % Execute D-ADMM
        
        % Optional input
        ops = struct('rho', {step_size}, ...
            'max_iter', {max_iter}, ...
            'x_opt', {Tx_psd'}, ...
            'eps_opt', {1e-2}, ...
            'turn_off_eps', {0}....
            );
        
        cd(direct_DADMM);
        tic
        [W, Y, vars_prob, ops_out] = DADMM(n, vars_prob, vars_network, ops);
        toc
        cd(direct_current);
        mse_soln(j) = norm(solution' - Tx_psd)/norm(Tx_psd);
        mse_admm(j) = norm(Z{1}' - Tx_psd)/norm(Tx_psd);
    end
    mse_admm_av(i) = sum(mse_admm)/runs;
    mse_soln_av(i) = sum(mse_soln)/runs;
end

iterations = ops_out.iterations;
error(1,:) = ops_out.error_iterations_z;
error_accel(1,:) = ops_out_accel.error_iterations_z;

figure
semilogy(1:iterations, error(1:iterations), 'r' , 1:iterations, error_accel(1:iterations), 'b');
legend('normal','accel')
