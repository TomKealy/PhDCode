clear all;
close all;

%randn('seed', 0);
%rand('seed',0);

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

Inc = adj2inc(Adj);
w_local_degree = local_degree(Adj);
%weights = construct_weight_mtx(w_local_degree, Adj);
%Create struct with network information
vars_network = struct('P', {P}, ...
    'neighbors', {neighbors}, ...
    'W', {w_local_degree}, ...
    'partition_colors', {partition_colors} ...
    );

m = 50;

L = 300;
M = L;

edges = [50, 120, 170, 192, 220, 224, 256, 300] ;
levels = [100,  0 , 400, 0, 0, 200, 0, 0];
idxs = zeros(1, M)  ;
idxs(edges(1: end-1)+1) = 1 ;
g = levels(cumsum(idxs)+1) ;

edges = [50, 120, 170, 192, 220, 234, 256, 300] ;
levels = [100,  0 , 400, 0, 0, 0, 0, 0];
idxs = zeros(1, M)  ;
idxs(edges(1: end-1)+1) = 1 ;
gstar = levels(cumsum(idxs)+1) ;

positions = randi(L,[1,5]);%generate random spikes for signal

Tx_psd = zeros(1,L); %Tx PSD
Tx_psd(positions) = 1000;

S = randn(m,L);
sigma = 10^(-10/20);
eta = randn(1,m)/m;
noise_sum = sum(eta);

Tx_psd = g;

A_BP = randn(m, M);

b = A_BP*Tx_psd';%+ sigma*eta';

lambda = sqrt(2*log(size(Tx_psd,2)));

max_iter = 200;
error_i_accel = zeros(1, max_iter+1);
path = zeros(100, L);

F = LehmerMatrix(M);
[l, U] = lu(F);

max_eig = max(abs(eig((A_BP*l)'*(A_BP*l))));

rho = 1/max_eig;
rho = nthroot(rho, 3);

if mod(m,P) ~= 0
    error('m divided by P must be integer');
end
m_p = m/P; 

% Create struct with problem data used in 'minimize_quad_prog_plus_l1_BB'
vars_prob = struct('handler_GPSR', @GPSR_BB, ...
    'A_BPDN', {A_BP}, ...
    'b_BPDN', {b}, ...
    'm_p', {m_p}, ...
    'Template', {Tx_psd'}, ...
    'D0', {l}, ...
    'P', {P}, ...
    'beta', {lambda*max_eig}, ...
    'relax', {1}...
    );

% Optional input
ops = struct('rho', {rho}, ...
    'max_iter', {max_iter}, ...
    'x_opt', {Tx_psd'}, ...
    'eps_opt', {1e-2}, ...
    'turn_off_eps', {0} ...
    );

cd(direct_DADMM);
tic
[X, Z, D, vars_prob, ops_out_lo] = DADMM_dl(L, vars_prob, vars_network, ops);
toc
cd(direct_current);

%============================================================================
