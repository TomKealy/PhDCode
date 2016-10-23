function [X, vars_prob, varargout] = DLasso(n, vars_prob, vars_network,...
    varargin) 

% [X, vars_prob, varargout] = DLasso(n, vars_prob, vars_network,...
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
% in a simulated distributed scenario using the D-Lasso
% 
%  [1] J. Bazerque and G. Giannakis, "Distributed Spectrum Sensing for 
%      Cognitive Radio Networks by exploiting Sparsity," IEEE Trans. Sig. 
%      Proc., Vol. 58, No. 3, 2010
%
%  Although we name the algorithm D-Lasso, it first appeared under the name
%  CA-MoM (consensus averaging method of multipliers) in
%
%  [2] H. Zhu, G. Giannakis, and A. Cano, "Distributed In-Network Channel
%      Decoding," IEEE Trans. Sig. Proc., Vol. 57, No. 10, 2009
%
%
% The size of the variable x is n. It is only required each function
% fp(.) to be convex and each set Xp to be convex. D-Lasso is oblivious to 
% the format of either fp or Xp. It only requires the user to provide a 
% function that solves, for each p,
%
%           minimize    fp(x) + v'*x + c||x||^2          (2)
%           subject to  x in Xp
%
% where the vector v and the scalar c are inputs given by D-Lasso. The 
% input vars_prob, which is a struct, should contain the handler of the
% function solving (2) in the field 'handler', i.e., 
% vars_prob.handler = function_handler. The remaining fields of vars_prob 
% are not used by D-Lasso and can be used to store internal information for 
% the solver of (2). The header of the solver of (2) should be:
%
%           [xp, vars_prob] = function_name(p, v, c, X, vars_prob)
%
% where p is the node number, v and c are as in (2), and X is a cell array
% such that X{p} is the last solution returned by this function, in xp, for
% node p. (Having X as an input reduces the memory usage whenever the 
% solver of (2) requires warm-starts: note that X is an internal variable 
% of D-Lasso.) The output xp is the solution of (2), and vars_prob is also
% returned because it might have changed (that depends on the user's 
% implementation).
%
% The remaining inputs and outputs of D-Lasso are now explained.
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
%                     network. 
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
%                           criteria are turned off: D-Lasso terminates 
%                           only when it reaches the maximum number of 
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
%                  may relevant information, D-Lasso returns it to the 
%                  user.
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
optargs = nargin - 3;
if optargs > 1
    error('Number of input arguments exceeds the maximum value. Please type ''help DLasso''.');
end
if optargs == 1
    if ~isstruct(varargin{1})
        error('Optional input should be a struct. Please type ''help DLasso''.');
    end
end
if ~isstruct(vars_prob)
    error('vars_prob (2nd argument of DLasso) should be a struct. Please type ''help DLasso''.');
end
if ~isfield(vars_prob, 'handler')
    error('vars_prob (2nd argument of DLasso) should contain the field ''handler''. Please type ''help DLasso''.');
end
if ~isstruct(vars_network)
    error('vars_network (3rd argument of DLasso) should be a struct. Please type ''help DLasso''.');
end
if sum(isfield(vars_network, {'P', 'neighbors'})) ~= 2
    error('Fieldnames of the struct vars_network are not correct. Please type ''help DLasso''.');
end
nout = max(nargout,1)-2;
if nout > 1
    error('Number of outputs is greater than the maximum. Please type ''help DLasso''.'); 
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
neighbors = vars_network.neighbors;
% =========================================================================


% =========================================================================
% Initializations

Lambda = cell(P,1);              % Cell (one/node) of dual variables
Grad_Lambda = cell(P,1);         % Cell (one/node) with X{i}-X{j}
X = cell(P,1);                   % Cell (one/node) with current estimates
X_prev = cell(P,1);              % Cell (one/node) with previous estimates

for p = 1 : P
    Lambda{p} = zeros(n,1);
    Grad_Lambda{p} = zeros(n,1);
    X{p} = zeros(n,1);
    X_prev{p} = zeros(n,1);
end

user_solver = vars_prob.handler; % Function provided by the user


Stop = zeros(P,1);               % if Stop(p) == 1, node p has terminated
                                     
                                      
% Struct with the internal variables of D-Lasso. 
vars = struct( 'X', {X},...
    'X_prev', {X_prev},...
    'EPS', {EPS},...
    'TURN_OFF_EPS', {TURN_OFF_EPS},...
    'EXISTS_X_OPT', {EXISTS_X_OPT},...
    'x_opt', {0},...                    % To be filled if EXISTS_X_OPT == 1
    'error_fun', {0},...                % To be filled if EXISTS_X_OPT == 1
    'error_fun_handler', {0},...        % To be filled if EXISTS_X_OPT == 1
    'EPS_OPT', {EPS_OPT} ...
    );                                      


if EXISTS_X_OPT == 1
    if ERROR_FUN == 1
        error_fun_handler = opt_input.error_fun;
        vars.error_fun = 1;
        vars.error_fun_handler = error_fun_handler;
        [error_iterations, vars_prob] = error_fun_handler(X, x_opt, vars_prob);
    else
        error_iterations = norm(X{1}-x_opt)/norm(x_opt);
    end
    iter_for_errors = zeros(10,2);
    iter_for_errors(:,1) = 10.^(-(1:10))';    % The first column has the errors
    iter_for_errors(:,2) = Inf*ones(10,1);    % The second column the number of iterations
    
    vars.x_opt = x_opt;
end
% =========================================================================
                                        
                          
% =========================================================================
% Algorithm

for k = 1 : MAX_ITER    

    % Update the primal variables x
    X_prev = X;
    X_aux = X; 
    for p = 1 : P
        
        neighbs = neighbors{p};
        Dp = length(neighbs);
        
        % Determining        
        % u = sum_{j in N_p} x
        u = zeros(n,1);
        for j = 1 : Dp
            u = u + X{neighbs(j)};
        end
        
        c = Dp*rho;
        v = 2*Lambda{p} - rho*Dp*X{p} - rho*u;
        [X_aux{p}, vars_prob] = user_solver(p, v, c, X, vars_prob);         
    end
    X = X_aux;
    
    
    % Update the dual variables Lambda
    for p = 1 : P
 
        neighbs = neighbors{p};
        Dp = length(neighbs);
        
        % Determining
        % u = sum_{j in N_p} x
        u = zeros(n,1);
       
        for j = 1 : Dp             
            u = u + X{neighbs(j)};   
        end
        
        Grad_Lambda{p} = (Dp*X{p} - u);
        
        Lambda{p} = Lambda{p} + (rho/2)*Grad_Lambda{p};
    end
    
                    
    if EXISTS_X_OPT == 1
        if ERROR_FUN == 1
            [new_error, vars_prob] = error_fun_handler(X, x_opt, vars_prob);
        else
            new_error = norm(X{1} - x_opt)/norm(x_opt);
        end
        error_iterations = [error_iterations, new_error]; %#ok<AGROW>
        ind_nonfilled = (iter_for_errors(:,2) == Inf);
        ind_lesserror = (new_error < iter_for_errors(:,1));
        intersect = ind_nonfilled & ind_lesserror;
        iter_for_errors(intersect,2) = k;
    end
    
    
    vars.X = X;
    vars.X_prev = X_prev;
    % Stopping criterion
    Stop = DLasso_stopping_criterion(Stop, vars_network, vars, vars_prob, k);
    
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
        error_iterations_out = error_iterations;
        iter_for_errors_out  = iter_for_errors;
    else
        error_iterations_out = [];
        iter_for_errors_out = [];
    end
    
    if k == MAX_ITER
        stop_crit = 'MAX ITERATIONS REACHED';
    else
        stop_crit = 'EPS REACHED';
    end
    
    varargout{1} = struct('iterations', k,...
        'stop_crit', stop_crit,...
        'error_iterations', error_iterations_out,...
        'iter_for_errors', iter_for_errors_out);
end
% =========================================================================


end



function [Stop] = DLasso_stopping_criterion(Stop, vars_network, vars, vars_prob, k)


if k <= 2 || vars.TURN_OFF_EPS == 1
    return
end

% Network variables
P = vars_network.P;

% D-Lasso variables
X = vars.X;
X_prev = vars.X_prev;
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
            new_error = norm(X{1} - x_opt)/norm(x_opt);
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
