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


% =========================================================================

M = 256;
K = 200;

edges = [50, 120, 170, 192, 220, 224, 256] ;
levels = [0,  0 , 400, 0, 0, 0, 0];
idxs = zeros(1, M)  ;
idxs(edges(1: end-1)+1) = 1 ;
Tx_psd = levels(cumsum(idxs)+1)' ;
npsd = Tx_psd;
sigma = 1;
noise = normrnd(0, 1, [1,M]);
Tx_psd = Tx_psd + noise';

H = haarmtx(M); % ConstructHaarWaveletTransformationMatrix(M) ;

%Tx_psd = H(2,:)' - 0.5*H(3,:)' + 0.5*H(4,:)' - H(5,:)' + H(6,:)' - 0.5*H(7,:)' + H(8,:)' - 0.25*H(9,:)' - 1.1*H(3,:)' + 0.5*H(51,:)' - 0.25*H(128,:)' - 1.2*H(35,:)' + H(3,:)' - 0.1*H(201,:)';

% % %d i f f e r e n t i a t i o n matrix
n(1:M-1) = -1;
D = diag(n, -1) + speye(M) ;

% %K Random measurements%
S = rand(K, M);
A_BP = S*H;
%s = spline(1:M, Tx_psd);
%w = wavedec(Tx_psd,8,'haar');       
b = A_BP*Tx_psd;
eg = max(abs(eig((A_BP)'*(A_BP))));
rho = nthroot(1/eg,3);
lambda = sqrt(2*log(M));
relaxation_parameter = 1.0;

[z0_r, history] = lasso_admm_1(S, b, lambda, rho, relaxation_parameter);

y0 = z0_r;

figure
plot(1:M, H*Tx_psd, 'b', 1:M, y0, 'm')
legend('true', 'estimate')
title('Normal Lasso')

% partition = [4, 1, 2, 4, 3, 9, 1, 2, 21, 7, 53, 1, 148];
% 
% [x history] = group_lasso(A_BP, b, lambda, partition, rho, 1.0);
% 
% figure
% plot(1:M, Tx_psd, 'b', 1:M, x, 'm')
% legend('true', 'estimate')
% title('Group Lasso')

% =========================================================================


% =========================================================================

max_iter = 10000;

%Check if matrix partition is possible (all blocks with same size)
if mod(K,P) ~= 0
    error('m divided by P must be integer');
end
m_p = K/P;                         % Number of rows of A each node stores

% Create struct with problem data used in 'minimize_quad_prog_plus_l1_BB'
vars_prob = struct('handler', @BPDN_RP_Solver,...
    'handler_GPSR', @GPSR_BB, ...
    'A_BPDN', {S}, ...
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
ops = struct('rho', {rho}, ...
    'max_iter', {max_iter}, ...
    'x_opt', {npsd}, ...
    'eps_opt', {1e-2}, ...
    'turn_off_eps', {0}....
    );

cd(direct_DADMM);
tic
[X, Z, vars_prob, ops_out_accel] = DADMM(M, vars_prob, vars_network, ops);
toc
cd(direct_current);

y1 = Z{22};

figure
plot(1:M, H*Tx_psd, 'b', 1:M, y1, 'm')
legend('true', 'estimate')
