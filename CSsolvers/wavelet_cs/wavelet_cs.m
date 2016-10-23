%direct_current = pwd;                                      % Current
%direct_networks = '/home/tk12098/Documents/MATLAB/DADMM_paper_Code/GenerateData/Networks';           % Networks

% =========================================================================

% =========================================================================

% Selecting the network

%cd(direct_networks);
load Nets_50_nodes.mat;     % File with randomly generated networks
%cd(direct_current);

% Select the network number: 1 to 7
net_num = 4 ;

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
levels = [300,  0 , 0, 400, 0, 0, 0];
idxs = zeros(1, M)  ;
idxs(edges(1: end-1)+1) = 1 ;
Tx_psd = levels(cumsum(idxs)+1)' ;
npsd = Tx_psd;
sigma = 1/M:64/M:5;
mse = zeros(1, size(sigma,2));
mse_av = zeros(1,100);

for j=1:1
    for i =1:1
        i
        noise = sigma(i)*normrnd(0, 1, [1,M]);
        Tx_psd = Tx_psd + noise';
        
        H = ConstructHaarWaveletTransformationMatrix(M) ;
        
        %Tx_psd = H(2,:)' - 0.5*H(3,:)' + 0.5*H(4,:)' - H(5,:)' + H(6,:)' - 0.5*H(7,:)' + H(8,:)' - 0.25*H(9,:)' - 1.1*H(3,:)' + 0.5*H(51,:)' - 0.25*H(128,:)' - 1.2*H(35,:)' + H(3,:)' - 0.1*H(201,:)';
        
        % % %d i f f e r e n t i a t i o n matrix
        n(1:M-1) = -1;
        D = diag(n, -1) + speye(M) ;
        
        % %K Random measurements%
        S = rand(K, M);
        A_BP = S*H;
        s = spline(1:M, Tx_psd);
        % w = wavedec(Tx_psd, 4, 'db10');
        % A_BP = rand(size(w,1), size(w,1));
        b = A_BP*Tx_psd;
        eg = max(abs(eig((A_BP)'*(A_BP))));
        rho = nthroot(1/eg,3);
        lambda = sigma(i)*20*sqrt(2*log(M));
        relaxation_parameter = 1;
        
        %    [z0_r, history] = lasso_admm_1(S, b, lambda, rho, relaxation_parameter);
        %
        %     y0 = z0_r;
        %
        %     figure
        %     plot(1:M, Tx_psd, 'b', 1:M, H'*y0, 'm')
        %     legend('true', 'estimate')
        %     title('Normal Lasso')
        
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
            'beta', {sqrt(40*log(M))}, ...
            'relax', {1.0}...
            );
        % =========================================================================
        
        % =========================================================================
        % Execute D-ADMM
        
        % Optional input
        ops = struct('rho', {0.5}, ...
            'max_iter', {max_iter}, ...
            'x_opt', {H*Tx_psd}, ...
            'eps_opt', {1e-2}, ...
            'turn_off_eps', {0}....
            );
        
        tic
        [X, Z, vars_prob, ops_out_accel] = DADMM(M, vars_prob, vars_network, ops);
        toc
        
        mse(i) = ops_out_accel.error_iterations_z(max_iter);
        
        % y1 = Z{1};
        %
        % figure
        % plot(1:M, Tx_psd, 'b', 1:M, H'*y1, 'm')
        % legend('true', 'estimate')
        % title('Distributed')
        %
        % figure
        % plot(ops_out_accel.error_iterations_z)
        % title('Convergence Lasso')
    end
    mse_av(j) = sum(mse)/100;
end
save('mse.mat', 'mse');