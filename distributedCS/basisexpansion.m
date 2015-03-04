clear all;
direct_current = pwd;                                      % Current
direct_networks = '/home/tk12098/Documents/MATLAB/DADMM_paper_Code/GenerateData/Networks';           % Networks
% Compressed Sensing Data
direct_data = '/home/tk12098/Documents/MATLAB/DADMM_paper_Code/GenerateData/ProblemData/CompressedSensing';
direct_DADMM = '/home/tk12098/Documents/MATLAB/DADMM_paper_Code/D-ADMM';

freq = 0:260;

network = gengeonet(50, 0.25);

ns = 1681;
candidates_x = 0:0.025:1;
candidates_y = 0:0.025:1;
[x,y] = meshgrid(candidates_x, candidates_y);
candidates = zeros(ns,2);
candidates(:,1) = x(:);
candidates(:,2) = y(:)';

P = 50;

pos = network.Positions;
g_sr = zeros(P,ns);
active_tx = candidates(25, :);
g_sr_act = zeros(1,P);

for j = 1:P
    g_sr_act(j) = exp(-norm(pos(j,:) - active_tx));
end

g_sr(:,25) = g_sr_act'; %matrxof channel co-effs in freq band 3

nb = 8;

beta_bs = zeros(1, ns*8); %in this example beta = 1 only for base 3
beta_bs(121*2+25) = 1000; %the base spanned by pulses(3, :) is 1, and there's only 1 Tx. 121*2+25

%phi_r = kron(g_sr, psi)*kron(psi, beta_bs)';

%note: g_sr should be nb*P rows and ns*nb columns

%g_srf = zeros(nb*P,ns*nb); %matrix of channel coeffs AT EACH FREQ
%g_srf(50:99,242:362) = g_sr;
A = exprnd(0.25, [nb*P, ns*nb]);

phi_rf = A*beta_bs';

y = phi_rf'; %+ 0.001*randn(1,400); %observations taken by the Rx

%LASSO
solution = GPSR_BB(y', A, 1000, 'Verbose', 1, 'StopCriterion', 3, 'ToleranceA', 1e-4);

%BPDN
[s, r, g, info] = spg_bp(A, y');

m=400;

%ADMM Soln
if mod(m,P) ~= 0
    error('m divided by P must be integer');
end
m_p = m/P;                         % Number of rows of A each node stores

%========================================================================
% ADMM for BP
%========================================================================

%ADMM Soln
if mod(m,P) ~= 0
    error('m divided by P must be integer');
end
m_p = m/P;                         % Number of rows of A each node stores

%The function (that should be built by the user of D-ADMM) %
%minimize_quad_prog_plus_l1_BB solves the problem
%min {(1/P)||x||_1 + v'*x + c||x||^2 : Ap x = bp}.
%We use duality to solve that problem and we need to initialize the dual
%variable before using the function. There is one dual variable per node.

dual_var = cell(P,1);

for p = 1 : P                      % Initialization with zeros
    dual_var{p} = zeros(m_p,1);
end

%Create struct with problem data used in 'minimize_quad_prog_plus_l1_BB'
vars_prob = struct('handler', @minimize_quad_prog_plus_l1_BB,...
    'A_BP', {A}, ...
    'b_BP', {y'},...
    'dual_var', {dual_var} ...
    );

ops = struct('rho', {0.1}, ...%penalty param
    'max_iter', {2000}, ...
    'eps', {1e-5}, ...  
    'turn_off_eps', {0} ...
    );
 %'turn_off_eps', {0} ...
 %'x_opt', {solution},...
  %  'eps_opt', {1e-3}, ... 

vars_network = struct('P', {P}, ...
    'neighbors', {network.Neighbors}, ...
    'partition_colors', {network.Partition} ...
    );

cd(direct_DADMM);
[X, vars_prob, ops_out] = DADMM(ns*nb, vars_prob, vars_network, ops);
cd(direct_current);

% ========================================================================
% ADMM for BPDN
% ========================================================================
% 
% The function (that should be built by the user of D-ADMM) %
% minimize_quad_prog_plus_l1_BB solves the problem
% min {(1/P)||x||_1 + v'*x + c||x||^2 : Ap x = bp}.
% We use duality to solve that problem and we need to initialize the dual
% variable before using the function. There is one dual variable per node.
% 
% dual_var = cell(P,1);
% 
% for p = 1 : P                      % Initialization with zeros
%     dual_var{p} = zeros(m_p,1);
% end
% 
% Create struct with problem data used in 'minimize_quad_prog_plus_l1_BB'
% vars_prob = struct('handler', @minimize_quad_prog_plus_l1_BB,...
%     'A_BP', {g_srf}, ...
%     'b_BP', {y'},...
%     'dual_var', {dual_var} ...
%     );
% 
% ops = struct('rho', {1}, ...
%     'max_iter', {1000}, ...
%     'eps', {1e-4}, ...
%     'x_opt', {beta_bs'},...
%     'eps_opt', {1e-1}, ...
%     'turn_off_eps', {0} ...
%     );
% 
% vars_network = struct('P', {P}, ...
%     'neighbors', {network.Neighbors}, ...
%     'partition_colors', {network.Partition} ...
%     );
% 
% cd(direct_DADMM);
% [X, vars_prob, ops_out] = DADMM(ns*nb, vars_prob, vars_network, ops);
% cd(direct_current);
% 
% ===================================================================
% ADMM for LASSO
% ===================================================================
% 
% Check if matrix partition is possible (all blocks with same size)
% if mod(m,P) ~= 0
%     error('m divided by P must be integer');
% end
% m_p = m/P;                         % Number of rows of A each node stores
% 
% Create struct with problem data used in 'minimize_quad_prog_plus_l1_BB'
% vars_prob = struct('handler', @BPDN_RP_Solver,...
%     'handler_GPSR', @GPSR_BB, ...
%     'A_BPDN', {g_srf}, ...
%     'b_BPDN', {y'}, ...
%     'm_p', {m_p}, ...
%     'P', {P}, ...
%     'beta', {0.3} ... %regularisation param
%     );
% =========================================================================
% 
% =========================================================================
% Execute D-ADMM
% 
% Optional input
% ops = struct('rho', {0.05}, ...
%     'max_iter', {1000}, ...
%     'x_opt', {solution}, ...
%     'eps_opt', {1e-2} ...
%     );
% 
% vars_network = struct('P', {P}, ...
%     'neighbors', {network.Neighbors}, ...
%     'partition_colors', {network.Partition} ...
%     );
% 
% cd(direct_DADMM);
% [X, vars_prob, ops_out] = DADMM(ns*nb, vars_prob, vars_network, ops);
% cd(direct_current);

iterations = ops_out.iterations;
stop_crit = ops_out.stop_crit;
error_iterations = ops_out.error_iterations;
iter_for_errors = ops_out.iter_for_errors;

fprintf('norm(A*X{1} - b)/norm(b) = %f\n',norm(A*X{1} - y)/norm(y));
fprintf('||X{1}||_1 = %f\n',norm(X{1},1));
fprintf('||X{1} - solutionl||/||solution|| = %f\n',norm(X{1}-solution)/norm(solution));
fprintf('Number of iterations = %d\n', iterations);
fprintf('stop_crit = %s\n', stop_crit);
fprintf('iter_for_errors = \n');
num_rows = size(iter_for_errors, 1);
for i_g = 1 : num_rows
    fprintf('%E    %d\n', iter_for_errors(i_g,1), iter_for_errors(i_g,2));
end

figure(1);clf;
semilogy(1:iterations,error_iterations(1:iterations), 'b');
title('error\_{iterations}');

figure(2)
plot(solution)

figure(3)
plot(X{1})
