function [X, vars_prob, varargout] = DLassoCompos(n, vars_prob, vars_network,...
    varargin) 

% [X, vars_prob, varargout] = DLassoCompos(n, vars_prob, vars_network,...
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
% where fp(x) = gp(x) + hp(x), in a simulated distributed scenario using 
% the D-Lasso algorithm proposed in (Algorithm 3)
% 
%  [1] G. Mateos, J. Bazerque and G. Giannakis, "Distributed Sparse Linear 
%      Linear Regression," IEEE Trans. Sig. Proc., Vol. 58, No. 10, 2010
%
% The size of the variable x is n. It is only required each function
% gp(.) and hp(.) to be convex and each set Xp to be convex. The constraint
% x in Xp will be associated to the function gp. This version of D-Lasso
% requires the user to provide two functions. One that solves, for each p,
%
%
%           minimize    gp(x) + v'*x + c||x||^2          (2)
%           subject to  x in Xp
%
% where the vector v and the scalar c are given by D-Lasso, and another 
% function that solves, for each p,
%
%           minimize    hp(x) + v'*x + c||x||^2          (3)
%              x
%
% Note that (3) is unconstrained (although it can be constrained implicitly
% depending on how hp is defined). The input vars_prob, which is a struct, 
% should contain the handler of the function solving (2) in the field 
% 'handler_cons', and the handler of the function solving (3) in the field
% 'handler_uncons', i.e., 
%    
%     vars_prob.handler_cons = function_handler_of_(2)
%     vars_prob.handler_uncons = function_handler_of_(3) 
%
% The remaining fields of vars_prob are not used and can be used to store 
% internal information for the solvers of (2) and (3). The header for each
% solver should be:
%
%           [xp, vars_prob] = function_name(p, v, c, X, vars_prob)
%
% where p is the node number, v and c are as in (2)-(3), and X is a cell
% array such that X{p} is the last solution returned by this function. 
% (Having X as an input reduces the memory usage whenever the solvers of 
% (2) and (3) require warm-starts: note that X is an internal variable of 
% DLassoCompos.) The output xp is the solution of either (2) or (3), and 
% vars_prob is also returned because it might have changed (that depends on
% the user's implementation).
%
% The remaining inputs and outputs of D-LassoMateos are now explained.
%
% Inputs:
%     - n: is the size of the variable x.
% 
%     - vars_prob: is a struct that contains at least the fields 
%                  'handler_cons', and 'handler_uncons', which contain the
%                  handlers for the solvers of (2) and (3), respectively. 
%                  The rest of vars_prob can be used by the user to store
%                  information (so that it does not get erased when Matlab 
%                  leaves that function).
%                  
%
%     - vars_network: is a struct containing information about the
%                     network. 
%           . vars_network.P is the number of nodes P.                      
%           . vars_network.neighbors is a cell array of the size of the 
%             number of nodes P, and each entry i contains the neighbors of
%             node i.           
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
%             [error, vars_prob] = name_of_function(X, Y, x_opt, vars_prob)
%
%                        where X is a cell of size P where X{p} is the
%                        estimate of node p to problem (1), and Y is the 
%                        same object, but associated to the function h. 
%                        x_opt is the optimal solution provided by the 
%                        user, and vars_prob is the struct also provided by
%                        the user. This function is called at the end of 
%                        each iteration (all nodes have updated their 
%                        estimates)
%
%           . eps_opt: only valid in case x_opt exists. This turns off
%                      the eps stopping criteria (based on two consecutive
%                      iterations in node 1), and replaces it by
%                      max{||X{1} - x_opt||,||Y{1} - x_opt||}/||x_opt|| 
%                                                               <= eps_opt.
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
% NOTE:  error_iterations and iter_for_errors contain the errors relative 
%        to the variable x. If the user provides its own error function (in 
%        error_fun), then these outputs will be in respect to that error.
%
% NOTE2: when error_fun is not provided, this function gives as error 
%
%                  max{error_on_X, error_on_Y}
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
    error('Number of input arguments exceeds the maximum value. Please type ''help DLassoCompos''.');
end
if optargs == 1
    if ~isstruct(varargin{1})
        error('Optional input should be a struct. Please type ''help DLassoCompos''.');
    end
end
if ~isstruct(vars_prob)
    error('vars_prob (2nd argument of DLassoCompos) should be a struct. Please type ''help DLassoCompos''.');
end
if ~isfield(vars_prob, {'handler_cons', 'handler_uncons'})
    error('vars_prob (2nd argument of DLassoCompos) should contain the fields ''handler_cons'' and ''handler_uncons''. Please type ''help DLassoCompos''.');
end
if ~isstruct(vars_network)
    error('vars_network (3rd argument of DLassoCompos) should be a struct. Please type ''help DLassoCompos''.');
end
if sum(isfield(vars_network, {'P', 'neighbors'})) ~= 2
    error('Fieldnames of the struct vars_network are not correct. Please type ''help DLassoCompos''.');
end
nout = max(nargout,1)-2;
if nout > 1
    error('Number of outputs is greater than the maximum. Please type ''help DLassoCompos''.'); 
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

Mu = cell(P,1);                  % Cell (one/node) of dual variables
Grad_Mu = cell(P,1);             % Cell (one/node) with DpX{p}-sum_j X{j}
X = cell(P,1);                   % Cell (one/node) with current estimates
X_prev = cell(P,1);              % Cell (one/node) with previous estimates
Y = cell(P,1);                   % Cell (one/node) with current estimates
Y_prev = cell(P,1);              % Cell (one/node) with previous estimates
Gamma = cell(P,1);               % Cell (one/node) with dual variables 
                                 % associated to the xp = yp

for p = 1 : P
    Mu{p} = zeros(n,1);
    Grad_Mu{p} = zeros(n,1);
    X{p} = zeros(n,1);
    X_prev{p} = zeros(n,1);
    Y{p} = X{p};
    Y_prev{p} = Y{p}; 
    Gamma{p} = zeros(n,1);
end

user_solver_g = vars_prob.handler_cons;     % Function provided by the user
user_solver_h = vars_prob.handler_uncons;   % Function provided by the user


Stop = zeros(P,1);               % if Stop(p) == 1, node p has terminated
                                     
                                      
% Struct with the internal variables of D-Lasso. 
vars = struct( 'X', {X}, ...
    'X_prev', {X_prev}, ...
    'Y', {Y}, ...
    'Y_prev', {Y_prev}, ...
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
        [error_iterations, vars_prob] = error_fun_handler(X, Y, x_opt, vars_prob);                
    else
        error_iterations = max(norm(X{1}-x_opt)/norm(x_opt),norm(X{1}-x_opt)/norm(x_opt));
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

    % Update the primal variables x and y
    Y_prev = Y;
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
        
        c = rho*(2*Dp+1)/2;
        v = Mu{p} + Gamma{p} - rho*(Y{p} + Dp*X{p} + u);
        [X_aux{p}, vars_prob] = user_solver_g(p, v, c, X, vars_prob);
        
        % Now, update variable y
        v_y = -Gamma{p} - rho*X_aux{p};
        c_y = rho/2;
        [Y{p}, vars_prob] = user_solver_h(p, v_y, c_y, Y, vars_prob);        
    end
    X = X_aux;
    
    
    % Update the dual variables Mu and Gamma
    for p = 1 : P
 
        neighbs = neighbors{p};
        Dp = length(neighbs);
        
        % Determining
        % u = sum_{j in N_p} x
        u = zeros(n,1);
       
        for j = 1 : Dp             
            u = u + X{neighbs(j)};   
        end
        
        Grad_Mu{p} = (Dp*X{p} - u);
        
        Mu{p} = Mu{p} + rho*Grad_Mu{p};
        
        Gamma{p} = Gamma{p} + rho*(X{p} - Y{p});
    end
    
                    
    if EXISTS_X_OPT == 1
        if ERROR_FUN == 1
            [new_error, vars_prob] = error_fun_handler(X, Y, x_opt, vars_prob);
        else
            new_error = max(norm(X{1} - x_opt)/norm(x_opt), norm(Y{1} - x_opt)/norm(x_opt));
        end
        error_iterations = [error_iterations, new_error]; %#ok<AGROW>
        ind_nonfilled = (iter_for_errors(:,2) == Inf);
        ind_lesserror = (new_error < iter_for_errors(:,1));
        intersect = ind_nonfilled & ind_lesserror;
        iter_for_errors(intersect,2) = k;
    end
    
                                     %figure(1);clf;semilogy(error_iterations);drawnow
    vars.X = X;
    vars.X_prev = X_prev;
    vars.Y = Y;
    vars.Y_prev = Y_prev;
    % Stopping criterion
    Stop = DLassoCompos_stopping_criterion(Stop, vars_network, vars, vars_prob, k);
    
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



function [Stop] = DLassoCompos_stopping_criterion(Stop, vars_network, vars, vars_prob, k)


if k <= 2 || vars.TURN_OFF_EPS == 1
    return
end

% Network variables
P = vars_network.P;

% D-Lasso variables
X = vars.X;
X_prev = vars.X_prev;
Y = vars.Y;
Y_prev = vars.Y_prev;
EPS = vars.EPS;
EPS_OPT = vars.EPS_OPT;

EXISTS_X_OPT = vars.EXISTS_X_OPT;
ERROR_FUN = vars.error_fun;

if EXISTS_X_OPT == 1
    x_opt = vars.x_opt;
       
    if ~isempty(EPS_OPT)
        if ERROR_FUN == 1
            error_fun_handler = vars.error_fun_handler;
            [new_error, vars_prob] = error_fun_handler(X, Y, x_opt, vars_prob);
        else
            new_error = max(norm(X{1} - x_opt)/norm(x_opt),norm(Y{1} - x_opt)/norm(x_opt));
        end
        if new_error <= EPS_OPT
            Stop = ones(P,1);
        end
        
        return;
    end
end

for p = 1 : P
    if max(norm(X{p} - X_prev{p})/norm(X_prev{p}),norm(Y{p} - Y_prev{p})/norm(Y_prev{p})) <= EPS
        Stop(p) = 1;
    end
    
end

end
