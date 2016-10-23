clear all;
close all;
clear all;

direct_current = pwd;                                      % Current
direct_networks = '/home/tk12098/Documents/MATLAB/DADMMPaperCode/GenerateData/Networks';           % Networks
% Compressed Sensing Data
direct_data = '/home/tk12098/Documents/MATLAB/DADMM_paper_Code/GenerateData/ProblemData/CompressedSensing';
direct_DADMM = '/home/tk12098/Documents/MATLAB/DADMMPaperCode/DADMM';
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
Eb_N0_dB = [10];
runs = 10;
[r,c] = size(Eb_N0_dB);

mse_admm_av = zeros(1,c);
mse_soln_av = zeros(1,c);

mse_admm = zeros(1,runs);
mse_soln = zeros(1,runs);

S = randn(m,L);
h = exprnd(0.15);
H = diag(h);
A_BP = S;
sigma = 10^(-Eb_N0_dB\20);
eta = randn(1,m)/m;
noise_sum = sum(eta);
b = A_BP*Tx_psd' + sigma*eta';

lambda = sigma*sqrt(2*log(size(Tx_psd,2)));
lambdas = 0.01:0.01:1;
% For the groundtruth, we use the spgl1 solver

solution = spgl1(A_BP, b, 0, 0.0001, []);

fprintf('||A_BP*solution-b|| = %E\n', norm(A_BP*solution-b));
fprintf('norm(solution,1) = %E\n', norm(solution,1));

% Check if matrix partition is possible (all blocks with same size)
if mod(m,P) ~= 0
    error('m divided by P must be integer');
end
m_p = m/P;                         % Number of rows of A each node stores

max_iter = 10000;
error_i = zeros(1, max_iter+1);
error_i_accel = zeros(1, max_iter+1);
path = zeros(100, n);

% Create struct with problem data used in 'minimize_quad_prog_plus_l1_BB'
vars_prob = struct(    'handler_GPSR', @GPSR_BB, ...
    'A_BPDN', {A_BP}, ...
    'b_BPDN', {b}, ...
    'm_p', {m_p}, ...
    'P', {P}, ...
    'beta', {lambda}, ...
    'relax', {1}...
    );
% =========================================================================

% =========================================================================
% Execute D-ADMM

% Optional input
ops = struct('rho', {0.1}, ...
    'max_iter', {max_iter}, ...
    'x_opt', {solution}, ...
    'eps_opt', {1e-2}, ...
    'turn_off_eps', {0}....
    );

cd(direct_DADMM);
tic
[X, Z, vars_prob, ops_out] = DADMM_AMP(n, m/P, vars_prob, vars_network, ops);
toc
cd(direct_current);

% vars_prob = struct('handler', @BPDN_RP_Solver,...
%     'handler_GPSR', @GPSR_BB, ...
%     'A_BPDN', {A_BP}, ...
%     'b_BPDN', {b}, ...
%     'm_p', {m_p}, ...
%     'P', {P}, ...
%     'beta', {0.1}, ...
%     'relax', {0.6}...
%     );
% 
% cd(direct_DADMM);
% tic
% [W, Y, vars_prob, ops_out_accel] = DADMM_lo(n, vars_prob, vars_network, ops);
% toc
% cd(direct_current);
% % =========================================================================

% =========================================================================
% Print results

error_i_accel = ops_out_accel.error_iterations_z;

iterations = ops_out.iterations;
stop_crit = ops_out_accel.stop_crit;
error_i = ops_out.error_iterations_z;

figure;clf;
semilogy(1:iterations,error_i(1:iterations), 'b', 1:iterations,error_i_accel(1,1:iterations), 'r');
legend('ell_1','ell_0');
title('error\_{iterations_z}');


figure;clf;
plot(1:iterations,error_i(1:iterations), 'b', 1:iterations,error_i_accel(1,1:iterations), 'r');
legend('ell_1','ell_0');
title('error\_{iterations_z}');


