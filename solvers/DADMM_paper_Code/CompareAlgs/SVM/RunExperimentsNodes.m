function [] = RunExperimentsNodes(number_of_nodes)

P = number_of_nodes;

% We will select the rho that is best for 1e-4 of accuracy
rhos = [1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2];

% Filename for saving data
FILE_SAVE = ['Results/SVM_', num2str(P), '_nodes_Iris.mat'];

% Filename for the networks
FILE_NETS = ['../../GenerateData/Networks/Nets_', num2str(P), '_nodes.mat'];

% Filename for the data
FILE_DATA = '../../GenerateData/ProblemData/SVMs/Iris_overlap.mat';

% =========================================================================
% Directories where files are located

% Algorithms to be compared
dir_DADMM        = '../../D-ADMM';
dir_DLASSO       = '../../D-Lasso';
dir_DMLE         = '../../D-MLE';

% Directory where the solvers are
dir_handlers = '../../D-ADMM/SVM';

% Current directory
dir_Curr = pwd;
% =========================================================================


load(FILE_NETS);               % Load file with several networks

% =========================================================================
% Execution parameters of the algorithms: 
% Stop when true solution is reached within eps_opt, or when 
% max_communications is reached
%
max_communications = 1e4;
eps_opt = 1e-3;           % Has to be one of {1e-1, 1e-2, ..., 1e-9, 1e-10}
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

beta = 1;

load(FILE_DATA);               % Load datafile: has solution and theta

X = A_mat;
Y_vec = D_vec;
[K , n] = size(A_mat);

% Check if all nodes can store the same ammount of data
if mod(K,P) ~= 0
    error('K divided by P must be integer');
end
m = K/P;                         % Number of rows of A each node stores
% =====================================================================


% =====================================================================
% Compute the solution with quadprog
YX = diag(Y_vec)*X;
Y1 = Y_vec;
Id_K = eye(K);

A_qp = [-YX , Y1 , -Id_K ; zeros(K,n+1) , -Id_K];
b_qp = [-ones(K,1) ; zeros(K,1)];
H_qp = zeros(n+1+K);
H_qp(1:n,1:n) = eye(n);
f_qp = [zeros(n+1,1) ; beta*ones(K,1)];

x_qp = quadprog(H_qp,f_qp,A_qp,b_qp);

solution = x_qp(1:n+1);
% =====================================================================


% =====================================================================
% Store data in structure

A_qp = cell(P,1);
b_qp = [-ones(m,1) ; zeros(m,1)];  % b_p is the same for all nodes
Xi_p = cell(P,1);                  % Holds the slack variables of each node

Id_m = eye(m);
for p = 1 : P
    ind = (1 + (p-1)*m) : p*m;
    YX_p = YX(ind,:);
    Y1_p = Y1(ind,:);
        
    A_qp{p} = [-YX_p , Y1_p , -Id_m ; zeros(m,n+1) , -Id_m];
    Xi_p{p} = zeros(m,1);
end

% Create struct with problem data
cd(dir_handlers);
vars_prob = struct('handler', @SVM_solver, ...
    'A_qp', {A_qp}, ...                             
    'b_qp', {b_qp}, ...                             
    'Xi_p', {Xi_p}, ...                             
    'beta', {beta}, ...                             
    'm', {m}, ...
    'n', {n}, ...
    'P', {P} ...    
    );
cd(dir_Curr);
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
        [X_DADMM, vars_prob_dummy, ops_out] = DADMM(n+1, vars_prob, vars_network, ops);
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
        [X_DLasso, vars_prob_dummy, ops_out] = DLasso(n+1, vars_prob, vars_network, ops);
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
        [X_DMLE, vars_prob_dummy, ops_out] = DMLE(n+1, vars_prob, vars_network, ops);
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
    
                
       % =================================================================
    % Save data
    
    save(FILE_SAVE, 'rhos', 'Num_Networks', 'i_num_net', 'n', 'm', ...
        'max_communications', 'eps_opt', ...
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