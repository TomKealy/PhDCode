function [X, Z, vars_prob, varargout] = DADMM_group_lasso(n, vars_prob, vars_network,...
    varargin) 

% [X, vars_prob, varargout] = DADMM(n, vars_prob, vars_network,...
%  varargin)
%
% Solves the optimization problem
%
%           minimize    f1(x) + f2(x) + ... + fP(x)      (1)
%           subject to  x in X1
%                       x in X2
%                        (...)
%                       x in XP
%
% in a simulated distributed scenario using the distributed alternating 
% direction method of multipliers: D-ADMM. D-ADMM solves (1) in a network
% of P nodes, where node p only knows the function fp and the set Xp.
% Detailed information about D-ADMM can be found in
%
%   J. Mota, J. Xavier, P. Aguiar, and M. Püschel, "D-ADMM: A 
%   Communication-Efficient Distributed Algorithm For Separable 
%   Optimization," Arxiv: http://arxiv.org/abs/1202.2805, 2013
%   To appear in IEEE Transactions on Signal Processing
%
% and
%
%   J. Mota, J. Xavier, P. Aguiar, and M. Püschel, "Distributed Basis 
%   Pursuit," IEEE Transactions Signal Processing, Vol. 60, Issue: 4, 2012
%   (Arxiv link: http://arxiv.org/abs/1009.1128)
% 
% In this documentation, the size of the variable x will be n. Each 
% function fp is assumed convex and each set Xp is also assumed convex. 
% D-ADMM requires the user to provide a function that solves, for each p,
%
%           minimize    fp(x) + v'*x + c||x||^2          (2)
%           subject to  x in Xp
%
% where the vector v and the scalar c are inputs given by D-ADMM. The 
% solver for (2) should be provided by the user through the input struct 
% 'vars_prob.' More concretely, vars_prob should be contain a field 
% 'handler' with a function handler for the solver of (2). That is, 
% vars_prob.handler = function_handler, where function_handler is a handler
% to a solver of (2), written by the user. The remaining fields of 
% vars_prob are not used by D-ADMM; rather, they can be used to store 
% internal information for the solver of (2). The header of the solver of 
% (2) should be:
%
%           [xp, vars_prob] = function_name(p, v, c, X, vars_prob)
%
% where p is the node number, v and c are as in (2), and X is a cell array
% such that X{p} is the last solution returned by this function, in xp, for
% node p. (Having X as an input reduces memory usage whenever the solver of
% (2) requires warm-starts: note that X is an internal variable 
% of D-ADMM.) The output xp is the solution of (2), and vars_prob is also
% returned because it might have changed (that depends on the user's 
% implementation).
%
% The remaining inputs and outputs of D-ADMM are now explained.
%
% Inputs:
%     - n: is the size of the variable x.
% 
%     - vars_prob: is a struct that contains at least the field 'handler',
%                  which contains the handler of user provided function 
%                  solving (2), i.e., vars_prob.handler = @function_name. 
%                  The rest of vars_prob can be used by the user to store
%                  information (so that it does not get erased when
%                  Matlab leaves that function).
%                  
%
%     - vars_network: is a struct containing information about the
%                     network: 
%           . vars_network.P is the number of nodes P.                      
%           . vars_network.neighbors is a cell array of the size of the 
%             number of nodes P, and each entry i contains the neighbors of
%             node i.           
%           . vars_network.partition_colors is a cell array of the size of
%             the number of (proper) colors of the network graph. Each 
%             entry contains a vector with the nodes that have that color.
%             For example, for a network with 6 nodes, 
%             {[1,2], [3,4,5], 6} would mean that we have three colors: 
%             nodes 1 and 2 get the first color, 3, 4, and 5 the second, 
%             and node 6 gets the third color.
%
%
%     - varargin: the optional argument is a struct which contains the
%                 following fields.
%           . rho:  a positive constant associated to the augmented 
%                   Lagrangian of f(x). The default value is 1.
%           . max_iter: maximum number of iterations. The default is 100.
%           . eps: the algorithm stops either because of max_iter or 
%                  because  ||X{1} - X_prev{1}||/||X_prev{1}|| <= eps , 
%                  i.e., if the relative difference between two consecutive
%                  iterations of node 1 is less or equal to eps. The 
%                  default value is 1e-4.
%           . x_opt: the optimal solution of (1). If the optimal solution 
%                    is passed, several (optional) outputs are activated.
%           . error_fun: only valid if x_opt exists. It is a function
%                        handler to assess the error, e.g., if the user 
%                        wants to see the error in the dual variable. The
%                        header should be:
%
%             [error, vars_prob] = name_of_function(X, x_opt, vars_prob)
%
%                        where X is a cell of size P where X{p} is the
%                        estimate of node p to problem (1), x_opt is the
%                        optimal solution provided by the user, and
%                        vars_prob is the struct also provided by the user.
%                        This function is called at the end of each
%                        iteration (all nodes have updated their estimates)
%
%           . eps_opt: only valid in case x_opt exists. This turns off
%                      the eps stopping criteria (based on two consecutive
%                      iterations in node 1), and replaces it by
%                      ||X{1} - x_opt||/||x_opt|| <= eps_opt.
%                      eps_opt has to be > 0.
%           . turn_off_eps: if this field is present (no matter what
%                           value), the eps or the eps_opt stopping 
%                           criteria are turned off: D-ADMM terminates only
%                           when it reaches the maximum number of 
%                           iterations.
%           
%
%
% Outputs:
%
%     - X: is a P x 1 cell array where the ith entry contains the solution
%          x, solving the optimization problem, for node i.
%
%     - vars_prob: is the struct that the user provides as input. Since it
%                  may relevant information, D-ADMM returns it to the user.
%
%     - varargout: (optional) it is a struct with the following fields.
%           . iterations: total number of iterations
%           . stop_crit: string with 'MAX ITERATIONS REACHED', if the 
%                        maximum number of iterations was reached, or
%                        'EPS REACHED', otherwise.
%           . error_iterations: this output is activated when x_opt is 
%                               passed as input. It is a vector with the
%                               size of the number of iterations, and
%                               contains the relative error in each entry 
%                               (corresponding to each iteration)
%           . iter_for_errors: activated when x_opt is passed. It is a
%                              matrix 2x10:
%                                          [  1e-1  ,  iter1 ]
%                                          [  1e-2  ,  iter2 ]
%                                          [       ...       ]
%                                          [  1e-10 , iter10 ]
%                              that contains the number of iterations to
%                              achieve a "canonical" relative error of
%                              1e-1, 1e-2, ..., 1e-10.
%
%
% *************************************************************************
% NOTE: error_iterations and iter_for_errors contain the errors relative to
%       the variable x. If the user provides its own error function (in 
%       error_fun), then these outputs will be in respect to that error.
% *************************************************************************
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Please send any questions, comments, or bug reports to joaomota@cmu.edu.
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


OPT_OUTPUT = 0;    % 0 if there is no optional output; 1 otherwise

% =========================================================================
% Check for input and output errors
optargs = nargin - 4;
if optargs > 1
    error('Number of input arguments exceeds the maximum value. Please type ''help DADMM''.');
end
if optargs == 1
    if ~isstruct(varargin{1})
        error('Optional input should be a struct. Please type ''help DADMM''.');
    end
end
if ~isstruct(vars_prob)
    error('vars_prob (2nd argument of DADMM) should be a struct. Please type ''help DADMM''.');
end
if ~isfield(vars_prob, 'handler')
    error('vars_prob (2nd argument of DADMM) should contain the field ''handler''. Please type ''help DADMM''.');
end
if ~isstruct(vars_network)
    error('vars_network (3rd argument of DADMM) should be a struct. Please type ''help DADMM''.');
end
if sum(isfield(vars_network, {'P', 'neighbors', 'partition_colors'})) ~= 3
    error('Fieldnames of the struct vars_network are not correct. Please type ''help DADMM''.');
end
if ~iscell(vars_network.partition_colors)
    error('vars_network.partition_colors is not a cell. Please type ''help DADMM''.');
end
nout = max(nargout,1)-3;
if nout > 1
    error('Number of outputs is greater than the maximum. Please type ''help DADMM''.'); 
end
if nout == 1
    OPT_OUTPUT = 1; 
end
% =========================================================================


% =========================================================================
% Take care of optional input

% Predefined and default variables
rho = 1;
MAX_ITER = 100;
EPS = 1e-4;

% Optional input
EXISTS_X_OPT = 0;     % 0 if there is no x_opt; 1 otherwise
TURN_OFF_EPS = 0;     % 1 if stopping criterion is maximum number of 
                      % iterations only; 0 otherwise
                      
EPS_OPT = 0;         % Equal to 0 if turned off. Oth., it has the eps_opt value 


if ~isempty(varargin)
    opt_input = varargin{1};
    
    if isfield(opt_input, 'rho')
        rho = opt_input.rho;
    end
    if isfield(opt_input, 'max_iter')
        MAX_ITER = opt_input.max_iter;
    end
    if isfield(opt_input, 'eps')
        EPS = opt_input.eps;
    end
    if isfield(opt_input, 'x_opt')
        x_opt = opt_input.x_opt;
        EXISTS_X_OPT = 1;
    end
    
    ERROR_FUN = 0;        % 1 if user provides error function
    
    if EXISTS_X_OPT == 1
        if isfield(opt_input, 'error_fun')
            ERROR_FUN = 1;
        end
        if isfield(opt_input, 'eps_opt')
            EPS_OPT = opt_input.eps_opt;
        end
    end
    if isfield(opt_input, 'turn_off_eps')
        TURN_OFF_EPS = 1;
    end
end

% Network variables
P = vars_network.P;

% =========================================================================


% =========================================================================
% Initializations

U = cell(P,1);  
U_aux = cell(P,1);

T = cell(P,1);           
T_aux = cell(P,1);                % Cell (one/node) with the dual variable
Diff = cell(P,1);                % Cell (one/node) with 
                                 %    Dp*X{p} - sum_{j in neighbs}X{j} 
X = cell(P,1);                   % Cell (one/node) with current estimates
X_prev = cell(P,1);              % Cell (one/node) with previous estimates

Z = cell(P,1);
Z_prev = cell(P,1);

L = cell(P,1);

Ap = cell(P,1);
bp = cell(P,1);
Atb = cell(P,1);
gamma0 = cell(P,1);
gamma1 = cell(P,1);

for p = 1 : P
    X{p} = zeros(n,1);
    X_prev{p} = X{p};
    Z{p} = zeros(n,1);
    Z_prev{p} = Z{p};
    U{p} = zeros(n,1);
    Diff{p} = zeros(n,1);
    T{p} = zeros(n,1);
    T_aux{p} = zeros(n,1);
    U_aux{p} = zeros(n,1);
    L{p} = zeros(n,1);
    Ap{p} = vars_prob.A_BPDN(p,:);
    bp{p} = vars_prob.b_BPDN(p,:);
    Atb{p} = Ap{p}'*bp{p};
end

Stop = zeros(P,1);               % if Stop(p) == 1, node p has terminated
                                     
                                      
% Struct with the internal variables of D-ADMM. 
vars = struct( 'U', {U}, ...   
    'T', {T}, ...
    'Diff', {Diff}, ...
    'X', {X},...
    'Z', {Z},...
    'X_prev', {X_prev},...
    'Z_prev', {Z_prev},...
    'rho', {rho}, ...
    'EPS', {EPS},...
    'MAX_ITER', {MAX_ITER},...
    'Stop', {Stop},...
    'TURN_OFF_EPS', {TURN_OFF_EPS},...    
    'EPS_OPT', {EPS_OPT},...
    'EXISTS_X_OPT', {EXISTS_X_OPT},...
    'L', {L}, ....
    'Ap', {Ap}, ...
    'bp', {bp}, ...
    'Atb', {Atb}, ...
    'x_opt', {0},...                    % To be filled if EXISTS_X_OPT == 1
    'error_fun', {0},...                % To be filled if EXISTS_X_OPT == 1
    'error_fun_handler', {0},...        % To be filled if EXISTS_X_OPT == 1
    'error_iterations_z', {[]},...        % To be filled if EXISTS_X_OPT == 1
    'error_iterations_x', {[]},...        % To be filled if EXISTS_X_OPT == 1
    'iter_for_errors', {0}, ...         % To be filled if EXISTS_X_OPT == 1
    'iterations', {0},...                % To be filled if EXISTS_X_OPT == 1
    'gamma0', {gamma0},...
    'gamma1', {gamma1},...
    'T_aux', {T_aux},...
    'U_aux', {U_aux}...
    );                                      


if EXISTS_X_OPT == 1
    if ERROR_FUN == 1
        error_fun_handler = opt_input.error_fun;
        vars.error_fun = 1;
        vars.error_fun_handler = error_fun_handler;
        [error_iterations_z, vars_prob] = error_fun_handler(Z, x_opt, vars_prob);
    else
        size(x_opt)
        size(Z{1})
        error_iterations_z = norm(Z{1}-x_opt)/norm(x_opt);
        %error_iterations_x = norm(X{1}-x_opt)/norm(x_opt);
    end
    iter_for_errors = zeros(10,2);
    iter_for_errors(:,1) = 10.^(-(1:10))';    % The first column has the errors
    iter_for_errors(:,2) = Inf*ones(10,1);    % The second column the number of iterations
    iterations = 0;
    
    vars.x_opt = x_opt;
    vars.error_iterations_z = error_iterations_z;
    %vars.error_iterations_x = error_iterations_x;
    vars.iter_for_errors = iter_for_errors;
    vars.iterations = iterations;
end
% =========================================================================
                                        
                          
% =========================================================================
% Algorithm
                                
for k = 1 : MAX_ITER    
    k
    vars.k = k;
    vars.X_prev = X;
    vars.Z_prev = Z;
    [vars, vars_prob] = DADMM_compute_gradient_lagr(vars, vars_prob, vars_network);
    Diff = vars.Diff;
    X = vars.X;
    Z = vars.Z;
    
    for p = 1 : P
        % Only iterate if the nodes are still active
        if Stop(p) == 0  
            U{p} = U{p} + rho*(Diff{p});

            T{p} = T{p} + rho*(X{p} - Z{p});
            
        end
    end
    vars.U = U;
    vars.T = T;
                                              
    % Stopping criterion for the outer loop
    Stop = DADMM_stopping_criterion(Stop, vars_network, vars, vars_prob, k);
    
    vars.Stop = Stop;
    
    % If all nodes have converged, stop the algorithm
    if sum(Stop) == P
        break;
    end   
end

% =========================================================================


% =========================================================================
% Optional output
if OPT_OUTPUT == 1
    if EXISTS_X_OPT == 1
        error_iterations_out_x = vars.error_iterations_x;
        error_iterations_out_z = vars.error_iterations_z;
        iter_for_errors_out  = vars.iter_for_errors;
    else
        error_iterations_out_x = [];
        error_iterations_out_z = [];
        iter_for_errors_out = [];
    end
    
    if k == MAX_ITER
        stop_crit = 'MAX ITERATIONS REACHED';
    else
        stop_crit = 'EPS REACHED';
    end
    
    varargout{1} = struct('iterations', k,...
        'stop_crit', stop_crit,...
        'error_iterations_x', error_iterations_out_x,...
        'error_iterations_z', error_iterations_out_z,...
        'iter_for_errors', iter_for_errors_out);
end
% =========================================================================
                               
                                   
end

function [vars, vars_prob] = DADMM_compute_gradient_lagr(vars, vars_prob, vars_network)

% This function updates the primal variables according to the ADMM 
% algorithm. It uses the graph coloring to select which nodes to update
% sequentially.

% DADMM variables
U = vars.U;
T = vars.T;
Diff = vars.Diff;
X = vars.X;
Z = vars.Z;
rho = vars.rho;
Stop = vars.Stop;
EXISTS_X_OPT = vars.EXISTS_X_OPT;
ERROR_FUN = vars.error_fun;
alpha = vars_prob.relax;
partition = vars_prob.partition;

if EXISTS_X_OPT == 1
    x_opt = vars.x_opt;
    error_iterations_z = vars.error_iterations_z;
    error_iterations_x = vars.error_iterations_x;
    iter_for_errors = vars.iter_for_errors;
    iterations = vars.iterations;
end

% Network variables
P = vars_network.P;
neighbors = vars_network.neighbors;
partition_colors = vars_network.partition_colors;    

% =========================================================================
% Algorithm

n = length(X{1});

num_colors = length(partition_colors);

% cumulative partition
cum_part = cumsum(partition);


for color = 1 : num_colors
    
    X_aux = X;
    Z_aux = Z;
    for p = partition_colors{color}
        
        if ~Stop(p)              % Only compute if not has not terminated
            neighbs = neighbors{p};
            Dp = length(neighbs);
            
            % Determine the sum of the X's of the neighbors       
            sum_neighbs = zeros(n,1);                        
            for j = 1 : Dp
                sum_neighbs = sum_neighbs + X{neighbs(j)};
            end
            
            % Solving each problem in the alternating direction minimization
            v = U{p} - rho*sum_neighbs;
            t = T{p};
            Ap = vars.Ap{p}; %pth row of A
            Atb = vars.Atb{p};
                                  
            zp = Z_aux{p};
                      
            %x-update
            q = (Ap'*Ap+rho*(Dp+1)*speye(size(Ap,2)));
          
            xp = (Atb+rho*zp - t - v);
            
            xp = q\xp;
            
            xp = alpha*xp + (1-alpha)*zp;
            
            zold = zp;
            start_ind = 1;
            x_hat = alpha*xp + (1-alpha)*zold;
            for i = 1:length(partition),
                sel = start_ind:cum_part(i);
                zp(sel) = shrinkage(x_hat(sel) + (t(sel)/rho), (vars_prob.beta/P)/rho);
                start_ind = cum_part(i) + 1;
            end
            
            %zp = shrinkage(xp + (t/rho), (vars_prob.beta/P)/(rho));
                   
            X_aux{p} = xp;
            Z_aux{p} = zp;
        end
    end
    X = X_aux;
    Z = Z_aux;
    
end
% =========================================================================

% =========================================================================
% Output

for p = 1 : P
    neighbs = neighbors{p};
    Dp = length(neighbs);
    
    % Determine the sum of the X's of the neighbors
    sum_neighbs = zeros(n,1);
    for j = 1 : Dp
        sum_neighbs = sum_neighbs + (X{neighbs(j)});
    end
    
    Diff{p} = Dp*(X{p}) - sum_neighbs;
end
vars.Diff = Diff;
vars.X = X;
vars.Z = Z;

if EXISTS_X_OPT == 1
    if ERROR_FUN == 1
        error_fun_handler = vars.error_fun_handler;
        [new_error_z, vars_prob] = error_fun_handler(Z, x_opt, vars_prob);
    else
        new_error_x = norm(X{1} - x_opt)/norm(x_opt);
        new_error_z = norm(Z{1} - x_opt)/norm(x_opt);
    end
    error_iterations_z = [error_iterations_z, new_error_z];    %figure(1);clf;semilogy(error_iterations);drawnow;
    error_iterations_x = [error_iterations_x, new_error_x];  
    ind_nonfilled = (iter_for_errors(:,2) == Inf);
    ind_lesserror = (new_error_z < iter_for_errors(:,1));
    intersect = ind_nonfilled & ind_lesserror;
    iter_for_errors(intersect,2) = iterations + 1;
        
    vars.error_iterations_z = error_iterations_z;
    vars.error_iterations_x = error_iterations_x;
    vars.iter_for_errors = iter_for_errors;
    vars.iterations = iterations + 1;
end
% =========================================================================

end


function [Stop] = DADMM_stopping_criterion(Stop, vars_network, vars, vars_prob, k)


if k <= 2 || vars.TURN_OFF_EPS == 1
    return
end

% Network variables
P = vars_network.P;

% D-ADMM variables
X = vars.Z;
Z= vars.Z;
X_prev = vars.X_prev;
Z_prev = vars.Z_prev;
EPS = vars.EPS;
EPS_OPT = vars.EPS_OPT;

EXISTS_X_OPT = vars.EXISTS_X_OPT;
ERROR_FUN = vars.error_fun;

if EXISTS_X_OPT == 1
    x_opt = vars.x_opt;
    
    
    if ~isempty(EPS_OPT)
        if ERROR_FUN == 1
            error_fun_handler = vars.error_fun_handler;
            [new_error, vars_prob] = error_fun_handler(X, x_opt, vars_prob);
        else
            new_error = norm(Z{1} - x_opt)/norm(x_opt);
        end
        if new_error <= EPS_OPT
            Stop = ones(P,1);
        end
        
        return;
    end
end

for p = 1 : P
    if norm(X{p} - X_prev{p})/norm(X_prev{p}) <= EPS
        Stop(p) = 1;
    end
    
end

end

function z = shrinkage(x, kappa)
    %z = sign(x)*subplus(x-kappa);
    z = max( 0, x - kappa ) - max( 0, -x - kappa );
end

function p = objective(A, b, lambda, x, z)
    p = ( 1/2*sum((A*x - b).^2) + lambda*norm(z,1));
end