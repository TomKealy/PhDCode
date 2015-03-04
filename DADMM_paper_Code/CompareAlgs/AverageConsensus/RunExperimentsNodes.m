function [] = RunExperimentsNodes(number_of_nodes)

P = number_of_nodes;

% We will select the rho that is best for 1e-4 of accuracy
rhos = [1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2];

% Filename for saving data
FILE_SAVE = ['Results/Cons_', num2str(P), '_nodes.mat'];

% Filename for the networks
FILE_NETS = ['../../GenerateData/Networks/Nets_', num2str(P), '_nodes.mat'];

% Filename for the data
FILE_DATA = ['../../GenerateData/ProblemData/Consensus/Data_Av_Consensus_', num2str(P), '_nodes.mat'];

% =========================================================================
% Directories where files are located

% Algorithms to be compared
dir_DADMM    = '../../D-ADMM';
dir_DLASSO   = '../../D-Lasso';
dir_DMLE     = '../../D-MLE';
dir_Oreshkin = '../../SpecificProblemSolvers/ConsensusAlgorithms';

% Directory where functions associated to consensus are
dir_handlers = '../../D-ADMM/AverageConsensus';

% Current directory
dir_Curr = pwd;
% =========================================================================


load(FILE_NETS);               % Load file with several networks

n = 1;                         % Size of the variable (scalar in this case)

% =========================================================================
% Execution parameters of the algorithms: 
% Stop when true solution is reached within eps_opt, or when 
% max_communications is reached
%
max_communications = 1e3;
eps_opt = 1e-4;           % Has to be one of {1e-1, 1e-2, ..., 1e-9, 1e-10}
% =========================================================================


pos_iter_errors = -log10(eps_opt);  %  Corresponding row in iter_for_errors
len_rhos = length(rhos);

Num_Networks = length(Networks);    % Number of networks

errors_DADMM    = cell(Num_Networks,1);
errors_DLasso   = cell(Num_Networks,1);
errors_DMLE     = cell(Num_Networks,1);
errors_Oreshkin = cell(Num_Networks,1);

best_rhos_DADMM  = cell(Num_Networks,1);
best_rhos_DLasso = cell(Num_Networks,1);
best_rhos_DMLE   = cell(Num_Networks,1);

iter_for_errors_DADMM    = zeros(Num_Networks,1);
iter_for_errors_DLasso   = zeros(Num_Networks,1);
iter_for_errors_DMLE     = zeros(Num_Networks,1);
iter_for_errors_Oreshkin = zeros(Num_Networks,1);

% =====================================================================
% Extract Data

load(FILE_DATA);               % Load datafile: has solution and theta
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
    
    % Struct with data for solving each node's problem
    cd(dir_handlers);
    vars_prob = struct('handler', @AverageConsensusSubProb,...
        'theta', {theta} ...
        );
    cd(dir_Curr);
    
    % Struct with network data (for D-ADMM, D-Lasso, and D-MLE)
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
            'error_fun', @RelativeError_Av_Consensus ...
            );
        cd(dir_Curr);
        
        fprintf('DADMM: start\n');
        cd(dir_DADMM);
        [X_DADMM, vars_prob, ops_out] = DADMM(n, vars_prob, vars_network, ops);
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
        cd(dir_handlers);
        ops = struct('rho', {rhos(i_rhos)}, ...
            'max_iter', {max_communications}, ...
            'x_opt', {solution}, ...
            'eps_opt', {eps_opt}, ...
            'error_fun', @RelativeError_Av_Consensus ...
            );
        cd(dir_Curr);
        
        fprintf('DLasso: start\n');
        cd(dir_DLASSO);
        [X_DLasso, vars_prob, ops_out] = DLasso(n, vars_prob, vars_network, ops);
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
        cd(dir_handlers);
        ops = struct('rho', {rhos(i_rhos)}, ...
            'max_iter', {max_communications}, ...
            'x_opt', {solution}, ...
            'eps_opt', {eps_opt}, ...
            'error_fun', @RelativeError_Av_Consensus ...
            );
        cd(dir_Curr);
        
        fprintf('DMLE: start\n');
        cd(dir_DMLE);
        [X_DMLE, vars_prob, ops_out] = DMLE(n, vars_prob, vars_network, ops);
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
    % Oreshkin
    
    % ============================================
    % Construct Metropolis-Hastings Weight Matrix
    W = zeros(P,P);
    
    for p = 1 : P
        neighbs = neighbors{p};
        Dp_p = Dp(p);
        quant = 0;
        for j_Dp = 1 : Dp_p
            neighbor = neighbs(j_Dp);
            max_deg = max(Dp(p), Dp(neighbor));
            W(p,neighbor) = 1/(max_deg + 1);
            W(neighbor,p) = W(p,neighbor);
            quant = quant + 1/(max_deg + 1);
        end
        W(p,p) = 1 - quant;
    end
    % ============================================
    
    eigen_vals = eig(W);
    l2 = eigen_vals(P-1);
    
    par1 = -0.5;
    par2 = 0;
    par3 = 3/2;
    
    weight = -(1/(par2 + (par3-1)*l2)^2)*( (par3-1)*l2^2 + par2*l2 + 2*par1 ...
        + 2*sqrt( par1^2 + par1*l2*(par2 + (par3-1)*l2) ) );
    
    
    % Create struct with problem data
    vars_prob = struct('theta', {theta}, ...
        'par1', {par1}, ...
        'par2', {par2}, ...
        'par3', {par3}, ...
        'weight', {weight}, ...
        'W', {W} ...
        );
    
    
    % Optional input
    cd(dir_handlers);
    ops = struct('max_iter', {max_communications}, ...
        'x_opt', {solution}, ...
        'eps_opt', {eps_opt}, ...
        'error_fun', @RelativeError_Av_Consensus ...
        );
    cd(dir_Curr);
            
    fprintf('Oreshkin: start\n');
    cd(dir_Oreshkin);
    [X_ConsOr, vars_prob, ops_out] = Oreshkin(n, vars_prob, vars_network, ops);
    cd(dir_Curr);
    fprintf('Oreshkin: finish\n');
    
    iterations = ops_out.iterations;
    stop_crit = ops_out.stop_crit;
    error_iterations = ops_out.error_iterations;
    iter_for_errors = ops_out.iter_for_errors;
    
    fprintf('Number of iterations = %d\n', iterations);
    fprintf('stop_crit = %s\n', stop_crit);
    for i_g = 1 : 6
        fprintf('%E    %d\n', iter_for_errors(i_g,1), iter_for_errors(i_g,2));
    end
    
    errors_Oreshkin{i_num_net} = error_iterations;
    iter_for_errors_Oreshkin(i_num_net) = iter_for_errors(pos_iter_errors,2);
    %******************************************************************
    
    
    % =================================================================
    % Save data
       
   save(FILE_SAVE, 'rhos', 'n', 'Num_Networks', 'i_num_net', ...
        'max_communications', 'eps_opt', ...
        'errors_DADMM', 'best_rhos_DADMM', 'iter_for_errors_DADMM', ...
        'errors_DLasso', 'best_rhos_DLasso', 'iter_for_errors_DLasso', ...
        'errors_DMLE', 'best_rhos_DMLE', 'iter_for_errors_DMLE', ...
        'errors_Oreshkin', 'iter_for_errors_Oreshkin' ...
        );
    % =================================================================
                
end

figure(1);clf;
semilogy(iter_for_errors_DADMM(1:i_num_net), 'bo-');
hold on;
semilogy(iter_for_errors_DLasso(1:i_num_net), 'rs-');
semilogy(iter_for_errors_DMLE(1:i_num_net), 'go-');
semilogy(iter_for_errors_Oreshkin(1:i_num_net), 'ko-');
drawnow;
ylabel('Communication steps');
xlabel('Network number');
legend('D-ADMM', 'D-Lasso', 'D-MLE', 'Oreshkin');
ylim([1,1e4]);
