clear all;
close all;

%direct_current = pwd;                                      % Current
%direct_networks = '/home/tk12098/Documents/MATLAB/DADMMPaperCode/GenerateData/Networks';           % Networks
% Compressed Sensing Data
%direct_data = '/home/tk12098/Documents/MATLAB/DADMMPaperCode/GenerateData/ProblemData/CompressedSensing';
%direct_DADMM = '/home/tk12098/Documents/MATLAB/DADMMPaperCode/DADMM';


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

%Create struct with network information
vars_network = struct('P', {P}, ...
    'neighbors', {neighbors}, ...
    'partition_colors', {partition_colors} ...
    );

%=========================================================================

%Matrices
M = 300;
K = 200;

A = randn(K,M);

F = LehmerMatrix(M);
[L U] = lu(F);

%setup

%Generate Signal
edges = [50, 120, 170, 192, 220, 224, 256, 300] ;
levels = [50,  0 , 400, 0, 0, 200, 0, 0];
idxs = zeros(1, M)  ;
idxs(edges(1: end-1)+1) = 1 ;
g = levels(cumsum(idxs)+1) ;

%Generate Noise
sigma_squared = 1;

noise = sigma_squared*randn(1, M);

ge = g + noise;

he = cumsum(ge)';

%h = Fa;

snr = norm(g)/norm(noise);

y = A*ge';

eg = max(abs(eig((A)'*(A))));
rho = nthroot(1/eg,3);
lambdas = [1 3, 7, 9, 11, 15, 21, 28 ,40];
lambda = sqrt(sigma_squared*2*log(M));
relaxation_parameter = 1.0;
cent_path = zeros(size(lambdas,2),M);
dist_path = zeros(size(lambdas,2),M);
dist_errors = zeros(size(lambdas,2), 2501)

for i=1:1 %size(lambdas,2)
    [ahat, history] = lasso_admm_1(A*L', y, lambdas(i), rho, relaxation_parameter);
    cent_path(i,:) = ahat;
    
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
        'beta', {lambdas(i)}, ...
        'relax', {relaxation_parameter}...
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
    
    [X, Z, vars_prob, ops_out_accel] = DADMM(M, vars_prob, vars_network, ops);
    dist_path(i,:) = Z{1};
    dist_errors(i,:) = ops_out_accel.error_iterations_z;
end

save('continuation_cent_dist.mat', 'dist_path', 'cent_path', 'dist_errors')