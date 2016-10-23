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

% direct_current = pwd;                                      % Current
% direct_networks = 'GenerateData/Networks';           % Networks
% % Compressed Sensing Data
% direct_data = 'GenerateData/ProblemData/CompressedSensing';
% direct_DADMM = '/home/tk12098/Documents/MATLAB/DADMM_paper_Code';

clear all;
direct_current = pwd;                                      % Current
direct_networks = '/home/tk12098/Documents/MATLAB/DADMM_paper_Code/GenerateData/Networks';           % Networks
% Compressed Sensing Data
direct_data = '/home/tk12098/Documents/MATLAB/DADMM_paper_Code/GenerateData/ProblemData/CompressedSensing';
direct_DADMM = '/home/tk12098/Documents/MATLAB/DADMM_paper_Code/D-ADMM';  % D-ADMM
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

n=200;
m=50;
L=n;
positions = randi(L,[1,5]);%generate random spikes for signal

randn('seed', 0);
rand('seed',0);

Tx_psd = zeros(1,L); %Tx PSD
Tx_psd(positions) = 100;
Eb_N0_dB = [5:20];
runs = 10;
[r,c] = size(Eb_N0_dB);

mse_admm = zeros(1,runs);
mse_soln = zeros(1,runs);

signal_energy = norm(Tx_psd)^2;

for i=1:1
    i
    for j=1:1
        j
        S = randn(m,L);
               
        h = sqrt(1/2)*(randn(1,L)+1j*randn(1,L)); %TOK
        H = diag(abs(h)); %TOK

        A_BP = S;%*H;
        sigma = 10^(-Eb_N0_dB(i)/20)
        eta = sqrt(sigma)*(randn(1,m)/m);%+ i*randn(1,m));
        noise_energy = norm(eta)^2;
        SNR = signal_energy/noise_energy
        b = A_BP*Tx_psd'+ eta';
           
        % For the groundtruth, we use the spgl1 solver
        
        opts = spgSetParms('verbosity', 0); 
        solution = spgl1(A_BP, b, 0, 0.1, []);
                      
        % Check if matrix partition is possible (all blocks with same size)
        if mod(m,P) ~= 0
            error('m divided by P must be integer');
        end
        m_p = m/P;                         % Number of rows of A each node stores
        
        vars_prob = struct('handler', @BPDN_RP_Solver,...
            'handler_GPSR', @GPSR_BB, ...
            'A_BPDN', {A_BP}, ...
            'b_BPDN', {b}, ...
            'm_p', {m_p}, ...
            'P', {P}, ...
            'beta', {0.1}, ...
            'alpha', {1}...
            );
        % =========================================================================
        
        % =========================================================================
        % Execute D-ADMM
        
        % Optional input
        ops = struct('rho', {0.5}, ...
            'max_iter', {1000}, ...
            'x_opt', {Tx_psd'}, ...
            'eps_opt', {1e-2}, ...
            'turn_off_eps', {0}....
            );
        
        % =========================================================================
        % Execute D-ADMM
        
        cd(direct_DADMM);
        [X, Z, vars_prob, ops_out] = DADMM_small(n, vars_prob, vars_network, ops);
        cd(direct_current);
        % =========================================================================
        
        % =========================================================================
        % Print results
        
        iterations = ops_out.iterations;
        stop_crit = ops_out.stop_crit;
        error_iterations = ops_out.error_iterations_z;
        iter_for_errors = ops_out.iter_for_errors;
        
        fprintf('norm(A*x_estimate - b)/norm(b) = %f\n',norm(A_BP*Z{1} - b)/norm(b));
        fprintf('||x_estimate||_1 = %f\n',norm(Z{1},1));
        fprintf('||x_estimate - solutionl||/||solution|| = %f\n',norm(Z{1} - solution)/norm(solution));
        fprintf('Number of iterations = %d\n', iterations);
        fprintf('stop_crit = %s\n', stop_crit);
        fprintf('iter_for_errors = \n');
        num_rows = size(iter_for_errors, 1);
        for i_g = 1 : num_rows
            fprintf('%E    %d\n', iter_for_errors(i_g,1), iter_for_errors(i_g,2));
        end
        
        mse_soln(j) = norm(solution' - Tx_psd)/norm(Tx_psd);
        mse_admm(j) = norm(Z{1}' - Tx_psd)/norm(Tx_psd);
        
        figure
        plot(Z{1})
        legend('Z')
        
        figure
        plot(solution)
        legend('spgl1')
    end
    mse_admm_av(i) = sum(mse_admm)/runs;
    mse_soln_av(i) = sum(mse_soln)/runs;
end

figure
plot(Eb_N0_dB, mse_soln_av, Eb_N0_dB, mse_admm_av)
legend('central', 'distributed');

