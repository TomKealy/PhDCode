function [] = RunExperimentsNodes(number_of_nodes)

P = number_of_nodes;

% We will select the rho that is best for 1e-4 of accuracy
rhos = [1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2];

% Filename for saving data
FILE_SAVE = ['Results/BPDN_RP_', num2str(P), '_nodes_902.mat'];

% Filename for the networks
FILE_NETS = ['../../GenerateData/Networks/Nets_', num2str(P), '_nodes.mat'];

% Filename for the data
FILE_DATA = '../../GenerateData/ProblemData/CompressedSensing/Id902_BPDN.mat';

% =========================================================================
% Directories where files are located

% Algorithms to be compared
dir_DADMM        = '../../D-ADMM';
dir_DLASSO       = '../../D-Lasso';
dir_DMLE         = '../../D-MLE';
dir_DLassoCompos = '../../SpecificProblemSolvers/D-Lasso-Compos';

% Directory where the solvers are
dir_handlers = '../../D-ADMM/BPDN-RowPartition';
dir_handlers_DLassoCompos = '../../SpecificProblemSolvers/D-Lasso-Compos/BPDN-RowPartition';

% Current directory
dir_Curr = pwd;
% =========================================================================


load(FILE_NETS);               % Load file with several networks

% =========================================================================
% Execution parameters of the algorithms: 
% Stop when true solution is reached within eps_opt, or when 
% max_communications is reached
%
max_communications = 2e3;
eps_opt = 1e-4;           % Has to be one of {1e-1, 1e-2, ..., 1e-9, 1e-10}
% =========================================================================


pos_iter_errors = -log10(eps_opt);  %  Corresponding row in iter_for_errors
len_rhos = length(rhos);

Num_Networks = length(Networks);    % Number of networks

errors_DADMM        = cell(Num_Networks,1);
errors_DLasso       = cell(Num_Networks,1);
errors_DMLE         = cell(Num_Networks,1);
errors_DLassoCompos = cell(Num_Networks,1);

best_rhos_DADMM        = cell(Num_Networks,1);
best_rhos_DLasso       = cell(Num_Networks,1);
best_rhos_DMLE         = cell(Num_Networks,1);
best_rhos_DLassoCompos = cell(Num_Networks,1);

iter_for_errors_DADMM        = zeros(Num_Networks,1);
iter_for_errors_DLasso       = zeros(Num_Networks,1);
iter_for_errors_DMLE         = zeros(Num_Networks,1);
iter_for_errors_DLassoCompos = zeros(Num_Networks,1);

% =====================================================================
% Extract Data

beta = 0.3;                    % Weight given to ||.||_1 in BPDN

load(FILE_DATA);               % Load datafile: has solution and theta

[m,n] = size(A_BP);
% =====================================================================


% =====================================================================
% Exact solution to BPDN

cd(dir_handlers);
[solution] = GPSR_BB(b, A_BP, beta, 'Verbose', 0, 'StopCriterion', 3, ...
    'ToleranceA', 1e-6);
cd(dir_Curr);

fprintf('||A_BP*solution-b|| = %E\n', norm(A_BP*solution-b));
fprintf('norm(solution,1) = %E\n', norm(solution,1));
% =====================================================================


for i_num_net = 1 : Num_Networks             
    
    % Gather network data
            
    Network = Networks{i_num_net};
    P = Network.P;
    Adj = Network.Adj;
    Type = Network.Type;
    Parameters = Network.Parameters;
    partition_colors = Network.Partition;
    Num_Colors = Network.Num_Colors;              
    neighbors = Network.Neighbors;                
    Dp = Network.Degrees;                         
    %OptimRho = Network.OptimRho_Cons;
    
    fprintf('\nNetwork type = %s\n', Type);            
            
    % =================================================================
    % Execute Algorithms
    
    if mod(m,P) ~= 0
        error('m divided by P must be integer');
    end
    m_p = m/P;

    % Struct with data for solving each node's problem
    cd(dir_handlers);
    vars_prob = struct('handler', @BPDN_RP_Solver,...
        'handler_GPSR', @GPSR_BB, ...
        'A_BPDN', {A_BP}, ...
        'b_BPDN', {b}, ...
        'm_p', {m_p}, ...
        'P', {P}, ...
        'beta', {beta} ...
        );
    cd(dir_Curr);
            
    % Struct with network data
    vars_network = struct('P', {P}, ...
        'neighbors', {neighbors}, ...
        'partition_colors', {partition_colors} ...
        );
    
    
    % *****************************************************************
    % D-ADMM            
    best_iter_DADMM = Inf;
    for i_rhos = 1 : len_rhos
        
        % Optional input
        ops = struct('rho', {rhos(i_rhos)}, ...
            'max_iter', {max_communications}, ...
            'x_opt', {solution}, ...
            'eps_opt', {eps_opt} ...
            );
        
        fprintf('DADMM: start\n');
        cd(dir_DADMM);
        [X_DADMM, vars_prob_dummy, ops_out] = DADMM(n, vars_prob, vars_network, ops);
        cd(dir_Curr);
        fprintf('DADMM: finish\n');
        
        iterations = ops_out.iterations;
        stop_crit = ops_out.stop_crit;
        error_iterations = ops_out.error_iterations;
        iter_for_errors = ops_out.iter_for_errors;
        
        fprintf('Number of iterations = %d\n', iterations);
        fprintf('stop_crit = %s\n', stop_crit);
        for i_g = 1 : 6
            fprintf('%E    %d\n', iter_for_errors(i_g,1), iter_for_errors(i_g,2));
        end
        
        if iter_for_errors(pos_iter_errors,2) < best_iter_DADMM
            best_iter_DADMM = iter_for_errors(pos_iter_errors,2);
            errors_DADMM{i_num_net} = error_iterations;
            best_rhos_DADMM{i_num_net} = rhos(i_rhos);
            iter_for_errors_DADMM(i_num_net) = iter_for_errors(pos_iter_errors,2);
        end
    end
    % *****************************************************************
    
    
    % *****************************************************************
    % D-Lasso
    best_iter_DLasso = Inf;
    for i_rhos = 1 : len_rhos
        
        % Optional input
        ops = struct('rho', {rhos(i_rhos)}, ...
            'max_iter', {max_communications}, ...
            'x_opt', {solution}, ...
            'eps_opt', {eps_opt} ...
            );
        
        fprintf('DLasso: start\n');
        cd(dir_DLASSO);
        [X_DLasso, vars_prob_dummy, ops_out] = DLasso(n, vars_prob, vars_network, ops);
        cd(dir_Curr);
        fprintf('DLasso: finish\n');
        
        iterations = ops_out.iterations;
        stop_crit = ops_out.stop_crit;
        error_iterations = ops_out.error_iterations;
        iter_for_errors = ops_out.iter_for_errors;
        
        fprintf('Number of iterations = %d\n', iterations);
        fprintf('stop_crit = %s\n', stop_crit);
        for i_g = 1 : 6
            fprintf('%E    %d\n', iter_for_errors(i_g,1), iter_for_errors(i_g,2));
        end
        
        if iter_for_errors(pos_iter_errors,2) < best_iter_DLasso
            best_iter_DLasso = iter_for_errors(pos_iter_errors,2);
            errors_DLasso{i_num_net} = error_iterations;
            best_rhos_DLasso{i_num_net} = rhos(i_rhos);
            iter_for_errors_DLasso(i_num_net) = iter_for_errors(pos_iter_errors,2);
        end
    end
    % *****************************************************************
    
    
    % *****************************************************************
    % DMLE
    best_iter_DMLE = Inf;
    for i_rhos = 1 : len_rhos
        
        % Optional input
        ops = struct('rho', {rhos(i_rhos)}, ...
            'max_iter', {max_communications}, ...
            'x_opt', {solution}, ...
            'eps_opt', {eps_opt} ...
            );
        
        fprintf('DMLE: start\n');
        cd(dir_DMLE);
        [X_DMLE, vars_prob_dummy, ops_out] = DMLE(n, vars_prob, vars_network, ops);
        cd(dir_Curr);
        fprintf('DMLE: finish\n');
        
        iterations = ops_out.iterations;
        stop_crit = ops_out.stop_crit;
        error_iterations = ops_out.error_iterations;
        iter_for_errors = ops_out.iter_for_errors;
        
        fprintf('Number of iterations = %d\n', iterations);
        fprintf('stop_crit = %s\n', stop_crit);
        for i_g = 1 : 6
            fprintf('%E    %d\n', iter_for_errors(i_g,1), iter_for_errors(i_g,2));
        end
        
        if iter_for_errors(pos_iter_errors,2) < best_iter_DMLE
            best_iter_DMLE = iter_for_errors(pos_iter_errors,2);
            errors_DMLE{i_num_net} = error_iterations;
            best_rhos_DMLE{i_num_net} = rhos(i_rhos);
            iter_for_errors_DMLE(i_num_net) = iter_for_errors(pos_iter_errors,2);
        end
    end
    % *****************************************************************
    
    %******************************************************************
    % D-Lasso-Compos
    best_iter_DLassoCompos = Inf;
    for i_rhos = 1 : len_rhos
        
        rho = rhos(i_rhos);
        
        bpAp = cell(P,1);       % Will contain Ap'*bp
        M_inv = cell(P,1);      % Will contain inv(Ap'*Ap + rho*eye(n))        
        for p = 1 : P
            Ap = A_BP(1+(p-1)*m_p : p*m_p , :);
            bp = b(1+(p-1)*m_p : p*m_p);
            bpAp{p} = Ap'*bp;
            M_inv{p} = inv(Ap'*Ap + rho*eye(n));
        end
        
        cd(dir_handlers_DLassoCompos);
        vars_prob_DLC = struct('handler_cons', @Solver_BPDN_g, ...
            'handler_uncons', @Solver_BPDN_h, ...
            'M_inv', {M_inv}, ...
            'Apbp', {bpAp}, ...
            'P', {P}, ...
            'beta', {beta} ...
            );
        cd(dir_Curr);
                
        % Optional input
        ops = struct('rho', {rho}, ...
            'max_iter', {max_communications}, ...
            'x_opt', {solution}, ...
            'eps_opt', {eps_opt} ...
            );
        
        fprintf('D-Lasso-Compos: start\n');
        cd(dir_DLassoCompos);
        [X_LassoCompos, vars_prob_DLC, ops_out] = DLassoCompos(n, vars_prob_DLC, vars_network, ops);
        cd(dir_Curr);
        fprintf('D-Lasso-Compos: finish\n');
        
        iterations = ops_out.iterations;
        stop_crit = ops_out.stop_crit;
        error_iterations = ops_out.error_iterations;
        iter_for_errors = ops_out.iter_for_errors;
        
        fprintf('Number of iterations = %d\n', iterations);
        fprintf('stop_crit = %s\n', stop_crit);
        for i_g = 1 : 6
            fprintf('%E    %d\n', iter_for_errors(i_g,1), iter_for_errors(i_g,2));
        end
        
        if iter_for_errors(pos_iter_errors,2) < best_iter_DLassoCompos
            best_iter_DLassoCompos = iter_for_errors(pos_iter_errors,2);
            errors_DLassoCompos{i_num_net} = error_iterations;
            best_rhos_DLassoCompos{i_num_net} = rhos(i_rhos);
            iter_for_errors_DLassoCompos(i_num_net) = iter_for_errors(pos_iter_errors,2);
        end
    end
    %******************************************************************
        
    
    % =================================================================
    % Save data
       
   save(FILE_SAVE, 'rhos', 'n', 'Num_Networks', 'i_num_net', ...
        'max_communications', 'eps_opt', ...
        'errors_DADMM', 'best_rhos_DADMM', 'iter_for_errors_DADMM', ...
        'errors_DLasso', 'best_rhos_DLasso', 'iter_for_errors_DLasso', ...
        'errors_DMLE', 'best_rhos_DMLE', 'iter_for_errors_DMLE', ...
        'errors_DLassoCompos', 'best_rhos_DLassoCompos', 'iter_for_errors_DLassoCompos' ...
        );
    % =================================================================
                
end

figure(1);clf;
semilogy(iter_for_errors_DADMM(1:i_num_net), 'bo-');
hold on;
semilogy(iter_for_errors_DLasso(1:i_num_net), 'rs-');
semilogy(iter_for_errors_DMLE(1:i_num_net), 'go-');
semilogy(iter_for_errors_DLassoCompos(1:i_num_net), 'ks-');
drawnow;
ylabel('Communication steps');
xlabel('Network number');
legend('D-ADMM', 'D-Lasso', 'D-MLE', 'DLassoCompos');
ylim([1,1e4]);

