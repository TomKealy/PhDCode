function [] = RunExperimentsNodes(number_of_nodes)

P = number_of_nodes;

% We will select the rho that is best for 1e-4 of accuracy
rhos = [1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2];

% Filename for saving data
FILE_SAVE = ['Results/LASSO_CP_', num2str(P), '_nodes_902.mat'];

% Filename for the networks
FILE_NETS = ['../../GenerateData/Networks/Nets_', num2str(P), '_nodes.mat'];

% Filename for the data
FILE_DATA = '../../GenerateData/ProblemData/CompressedSensing/Id902_LASSO.mat';

% =========================================================================
% Directories where files are located

% Algorithms to be compared
dir_DADMM        = '../../D-ADMM';
dir_DLASSO       = '../../D-Lasso';
dir_DMLE         = '../../D-MLE';

% Directory where the solvers are
dir_handlers = '../../D-ADMM/LASSO-CP';

% Current directory
dir_Curr = pwd;
% =========================================================================


load(FILE_NETS);               % Load file with several networks

% =========================================================================
% Execution parameters of the algorithms: 
% Stop when true solution is reached within eps_opt, or when 
% max_communications is reached
%
max_communications = 1e3;
eps_opt = 5e-3;           % Has to be one of {1e-1, 1e-2, ..., 1e-9, 1e-10}
% =========================================================================


pos_iter_errors = -log10(eps_opt);  %  Corresponding row in iter_for_errors
len_rhos = length(rhos);

Num_Networks = length(Networks);    % Number of networks

errors_DADMM        = cell(Num_Networks,1);
errors_DLasso       = cell(Num_Networks,1);
errors_DMLE         = cell(Num_Networks,1);

best_rhos_DADMM        = cell(Num_Networks,1);
best_rhos_DLasso       = cell(Num_Networks,1);
best_rhos_DMLE         = cell(Num_Networks,1);

iter_for_errors_DADMM        = zeros(Num_Networks,1);
iter_for_errors_DLasso       = zeros(Num_Networks,1);
iter_for_errors_DMLE         = zeros(Num_Networks,1);

% =====================================================================
% Extract Data

delta = 1e-2;            % Regularization factor
sigma = 0.1;             % Parameter of LASSO

load(FILE_DATA);         % Load datafile

% Note that b was changed so that:
%b = 100*b + 0.1*sigma*randn(m,1);

[m,n] = size(A_BP);
% =====================================================================


% =====================================================================
% Exact solution to LASSO

options_spgl1 = struct('verbosity', {0});
solution = spgl1(A_BP, b, 0, sigma, [], options_spgl1);

fprintf('||A_BP*solution-b|| = %E < %E\n', norm(A_BP*solution-b), sigma);
fprintf('norm(solution,1) = %E\n', norm(solution,1));
fprintf('card(solution) = %d\n', sum(solution~=0));
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
    
    if mod(n,P) ~= 0
        error('n divided by P must be integer');
    end
    n_p = n/P;
    
    A = cell(P,1);
    max_sing_val_sq = cell(P,1);
    previous_t = cell(P,1);
    x_opt = cell(P,1);             % Estimates of x along the algorithm
    
    for p = 1 : P
        A{p} = A_BP(: , (p-1)*n_p + 1 : p*n_p );
        sing_vals = svd(A{p});
        max_sing_val_sq{p} = max(sing_vals)^2;
        previous_t{p} = 0;
        x_opt{p} = zeros(n_p,1);
    end
    
    one_over_P_b = 1/P*b;
    sig_P = sigma/P;
    
    
    % Create struct with problem data
    cd(dir_handlers);
    vars_prob = struct('handler', @NesterovForLassoSubProb,...
        'n', {n_p}, ...
        'A', {A}, ...
        'max_sing_val_sq', {max_sing_val_sq}, ...
        'delta', {delta}, ...
        'one_over_P_b', {one_over_P_b}, ...
        'sig_P', {sig_P}, ...
        'x_opt', {x_opt}, ...
        'previous_t', {previous_t} ...
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
        cd(dir_handlers);
        ops = struct('rho', {rhos(i_rhos)}, ...
            'max_iter', {max_communications}, ...
            'x_opt', {solution}, ...
            'eps_opt', {eps_opt}, ...
            'error_fun', @error_function_Lasso_for_DADMM ...
            );
        cd(dir_Curr);
        
        fprintf('DADMM: start\n');
        cd(dir_DADMM);
        [X_DADMM, vars_prob_dummy, ops_out] = DADMM(m, vars_prob, vars_network, ops);
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
        
        if iterations < best_iter_DADMM
            best_iter_DADMM = iterations;
            errors_DADMM{i_num_net} = error_iterations;
            best_rhos_DADMM{i_num_net} = rhos(i_rhos);
            iter_for_errors_DADMM(i_num_net) = iterations;
        end
    end
    % *****************************************************************
    
    
    % *****************************************************************
    % D-Lasso
    best_iter_DLasso = Inf;
    for i_rhos = 1 : len_rhos
        
        % Optional input
        cd(dir_handlers);
        ops = struct('rho', {rhos(i_rhos)}, ...
            'max_iter', {max_communications}, ...
            'x_opt', {solution}, ...
            'eps_opt', {eps_opt}, ...
            'error_fun', @error_function_Lasso_for_DADMM ...
            );
        cd(dir_Curr);
        
        fprintf('DLasso: start\n');
        cd(dir_DLASSO);
        [X_DLasso, vars_prob_dummy, ops_out] = DLasso(m, vars_prob, vars_network, ops);
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
        
        if iterations < best_iter_DLasso
            best_iter_DLasso = iterations;
            errors_DLasso{i_num_net} = error_iterations;
            best_rhos_DLasso{i_num_net} = rhos(i_rhos);
            iter_for_errors_DLasso(i_num_net) = iterations;
        end
    end
    % *****************************************************************
    
    
    % *****************************************************************
    % DMLE
    best_iter_DMLE = Inf;
    for i_rhos = 1 : len_rhos
        
        % Optional input
        cd(dir_handlers);
        ops = struct('rho', {rhos(i_rhos)}, ...
            'max_iter', {max_communications}, ...
            'x_opt', {solution}, ...
            'eps_opt', {eps_opt}, ...
            'error_fun', @error_function_Lasso_for_DADMM ...
            );
        cd(dir_Curr);
        
        fprintf('DMLE: start\n');
        cd(dir_DMLE);
        [X_DMLE, vars_prob_dummy, ops_out] = DMLE(m, vars_prob, vars_network, ops);
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
        
        if iterations < best_iter_DMLE
            best_iter_DMLE = iterations;
            errors_DMLE{i_num_net} = error_iterations;
            best_rhos_DMLE{i_num_net} = rhos(i_rhos);
            iter_for_errors_DMLE(i_num_net) = iterations;
        end
    end
    % *****************************************************************
    
    
        
    % =================================================================
    % Save data
       
   save(FILE_SAVE, 'rhos', 'Num_Networks', 'i_num_net', 'n', 'm', ...
        'max_communications', 'eps_opt', 'delta', 'sigma', ...
        'errors_DADMM', 'best_rhos_DADMM', 'iter_for_errors_DADMM', ...
        'errors_DLasso', 'best_rhos_DLasso', 'iter_for_errors_DLasso', ...
        'errors_DMLE', 'best_rhos_DMLE', 'iter_for_errors_DMLE' ...
        );
    % =================================================================
                
end

figure(1);clf;
semilogy(iter_for_errors_DADMM(1:i_num_net), 'bo-');
hold on;
semilogy(iter_for_errors_DLasso(1:i_num_net), 'rs-');
semilogy(iter_for_errors_DMLE(1:i_num_net), 'go-');
drawnow;
ylabel('Communication steps');
xlabel('Network number');
legend('D-ADMM', 'D-Lasso', 'D-MLE');
ylim([1,1e4]);
