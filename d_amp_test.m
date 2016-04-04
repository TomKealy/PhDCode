clear all;
close all;

randn('seed', 0);
rand('seed',0);


direct_current = pwd;                                      % Current
direct_networks = '/home/tk12098/Documents/MATLAB/DADMMPaperCode/GenerateData/Networks';           % Networks
% Compressed Sensing Data
%direct_data = '/home/tk12098/Documents/MATLAB/DADMMPaperCode/GenerateData/ProblemData/CompressedSensing';
direct_DADMM = '/home/tk12098/Documents/MATLAB/DADMMPaperCode/DADMM';


% =========================================================================


% =========================================================================

% Selecting the network

%cd(direct_networks);
load Nets_50_nodes.mat;     % File with randomly generated networks
%cd(direct_current);

% Select the network number: 1 to 7
net_num = 4 ;
%net_num = gengeonet(50, 0.25);

Adj = Networks{net_num}.Adj;                   % Adjacency matrix
partition_colors = Networks{net_num}.Partition;% Color partition of network

P = length(Adj);                               % Number of nodes

% Construct the cell neighbors, where neighbors{p} is a vector with the
% neighbors of node p
neighbors = cell(P,1);

for p = 1 : P
    neighbors{p} = find(Adj(p,:));
end

w_local_degree

%Create struct with network information
vars_network = struct('P', {P}, ...
    'neighbors', {neighbors}, ...
    'partition_colors', {partition_colors} ...
    );

L = 200;
m = 50;

positions = randi(L,[1,1]);%generate random spikes for signal

Tx_psd = zeros(1,L); %Tx PSD
Tx_psd(positions) = 1000;

S = randn(L,L);
A_BP = S/norm(S);
%sigma = 10^(-Eb_N0_dB\20);
eta = randn(1,m)/m;
noise_sum = sum(eta);
b = A_BP*Tx_psd';%$ + sigma*eta';

lambda = sqrt(2*log(size(Tx_psd,2)));

max_iter = 500;
error_i_accel = zeros(1, max_iter+1);
path = zeros(100, L);

max_eig = max(abs(eig(A_BP'*A_BP)));

rho = nthroot(1/max_eig, 3);

if mod(m,P) ~= 0
    error('m divided by P must be integer');
end
m_p = 4 ;%m/P; 

% Create struct with problem data used in 'minimize_quad_prog_plus_l1_BB'
vars_prob = struct('handler_GPSR', @GPSR_BB, ...
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
ops = struct('rho', {2.0}, ...
    'max_iter', {max_iter}, ...
    'x_opt', {Tx_psd'}, ...
    'eps_opt', {1e-2} ...
    );

cd(direct_DADMM);
tic
[X, Z, vars_prob, ops_out] = DADMM_AMP(L, vars_prob, vars_network, ops);
toc
cd(direct_current);

% % =========================================================================

% =========================================================================
% Print results

iterations = ops_out.iterations;
stop_crit = ops_out.stop_crit;
error_i = ops_out.error_iterations_x;

solution = spgl1(A_BP, b, 0, 0.0001, []);

fprintf('||A_BP*solution-b|| = %E\n', norm(A_BP*solution-b));
fprintf('norm(solution,1) = %E\n', norm(solution,1));

figure;clf;
semilogy(1:iterations,error_i(1:iterations), 'b');
legend('DAMP');
title('error\_{iterations_z}');

figure
plot(Z{1})

figure
plot(solution)
