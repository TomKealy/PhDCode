clear all;
close all;

direct_current = pwd;                                      % Current
direct_networks = '/home/tk12098/Documents/MATLAB/DADMMPaperCode/GenerateData/Networks';           % Networks
% Compressed Sensing Data
direct_data = '/home/tk12098/Documents/MATLAB/DADMMPaperCode/GenerateData/ProblemData/CompressedSensing';
direct_DADMM = '/home/tk12098/Documents/MATLAB/DADMMPaperCode/DADMM';


% =========================================================================


% =========================================================================

% Selecting the network

cd(direct_networks);
load Nets_50_nodes.mat;     % File with randomly generated networks
cd(direct_current);

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

% Create struct with network information
vars_network = struct('P', {P}, ...
    'neighbors', {neighbors}, ...
    'partition_colors', {partition_colors} ...
    );

% =========================================================================

data = csvread('tvws_data1.csv');
 
data_length = 8000;
ge = data(1:data_length,2);

M = length(ge);
K = 0.75*M;

F = LehmerMatrix(M);
[L U] = lu(F);

h = cumsum(g)';

he = cumsum(ge);

%h = Fa;

A = randn(K,M);

snr = norm(g)/norm(noise);

y = A*ge';

eg = max(abs(eig((A)'*(A))));
rho = nthroot(1/eg,3);
lambda = sqrt(12*2*log(M));
relaxation_parameter = 1.0;
path = zeros(size(lambdas,2),M);

%for i=1:size(lambdas,2)
    [ahat, history] = lasso_admm_1(A*L', y, lambda, rho, relaxation_parameter);
    %path(i,:) = ahat;
%end

str = sprintf('New basis, lambda = %f', lambda);

% figure
% plot(1:M, ge, 'b', 1:M, ahat, 'm')
% legend('orignal', 'denoised')
title(str)

max_iter = 2500;

%Check if matrix partition is possible (all blocks with same size)
if mod(K,P) ~= 0
    error('m divided by P must be integer');
end
m_p = K/P;                         % Number of rows of A each node stores

% Create struct with problem data used in 'minimize_quad_prog_plus_l1_BB'
vars_prob = struct('handler', @BPDN_RP_Solver,...
    'handler_GPSR', @GPSR_BB, ...
    'A_BPDN', {A*L'}, ...
    'b_BPDN', {y}, ...
    'm_p', {m_p}, ...
    'P', {P}, ...
    'beta', {lambda}, ...
    'relax', {1}...
    );
% =========================================================================

% =========================================================================
% Execute D-ADMM

% Optional input
ops = struct('rho', {rho}, ...
    'max_iter', {max_iter}, ...
    'x_opt', {ahat}, ...
    'eps_opt', {1e-2}, ...
    'turn_off_eps', {0}....
    );

cd(direct_DADMM);
tic
[X, Z, vars_prob, ops_out_accel] = DADMM(M, vars_prob, vars_network, ops);
toc
cd(direct_current);

y1 = Z{22};

save('new_basis_ofcom.mat', 'y1');

figure
plot(1:M, y1, 'b', 1:M, ahat, 'm')
legend('orignal', 'denoised')
title(str)