function [X, vars_prob, varargout] = Oreshkin(n, vars_prob, vars_network,...
    varargin) 

% [X, vars_prob, varargout] = Oreshkin(n, vars_prob, vars_network,...
%   varargin) 
%
% Implements the average consensus algorithm in 
%    B. Oreshkin, M. Coates, M. Rabbat,
%    "Optimization and Analysis of Distributed Averaging With Short Node
%    Memory," 
%    IEEE Trans. Sig. Proc., Vol. 58, No. 5, 2010.




OPT_OUTPUT = 0;    % 0 if there is no optional output; 1 otherwise

% =========================================================================
% Check for input and output errors
optargs = nargin - 3;
if optargs > 1
    error('Number of input arguments exceeds the maximum value. Please type ''help Oreshkin''.');
end
if optargs == 1
    if ~isstruct(varargin{1})
        error('Optional input should be a struct. Please type ''help Oreshkin''.');
    end
end
if ~isstruct(vars_prob)
    error('vars_prob (2nd argument of Oreshkin) should be a struct. Please type ''help Oreshkin''.');
end
if ~isstruct(vars_network)
    error('vars_network (3rd argument of Oreshkin) should be a struct. Please type ''help Oreshkin''.');
end
if sum(isfield(vars_network, {'P', 'neighbors'})) ~= 2
    error('Fieldnames of the struct vars_network are not correct. Please type ''help Oreshkin''.');
end
nout = max(nargout,1)-2;
if nout > 1
    error('Number of outputs is greater than the maximum. Please type ''help Oreshkin''.'); 
end
if nout == 1
    OPT_OUTPUT = 1; 
end
% =========================================================================


% =========================================================================
% Take care of optional input

% Predefined and default variables
MAX_ITER = 100;
EPS = 1e-4;

% Optional input
EXISTS_X_OPT = 0;     % 0 if there is no x_opt; 1 otherwise
TURN_OFF_EPS = 0;     % 1 if stopping criterion is maximum number of 
                      % iterations only; 0 otherwise
                      
EPS_OPT = 0;         % Equal to 0 if turned off. Oth., it has the eps_opt value 
                      
opt_input = varargin{1};

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

% Network variables
P = vars_network.P;
neighbors = vars_network.neighbors;
% =========================================================================


% =========================================================================
% Initializations

X = cell(P,1);                   % Cell (one/node) with current estimates
Z = cell(P,1);                   % Cell (one/node) with current estimates
Y = cell(P,1);                   % Cell (one/node) with current estimates
X_prev = cell(P,1);              % Cell (one/node) with previous estimates

theta = vars_prob.theta;         % Vector with initial estimates

par1 = vars_prob.par1;
par2 = vars_prob.par2;
par3 = vars_prob.par3;
weight = vars_prob.weight;
W = vars_prob.W;


Degrees = zeros(P,1);                 % Vector of degrees +1
for p = 1 : P
    neighbs = neighbors{p};
    Degrees(p) = length(neighbs) + 1;
end


for p = 1 : P
    Y{p} = zeros(n,1);
    Z{p} = zeros(n,1);
    X{p} = theta(p);
    X_prev{p} = X{p};
end

Stop = zeros(P,1);               % if Stop(p) == 1, node p has terminated
                                                                           

if EXISTS_X_OPT == 1
    if ERROR_FUN == 1
        error_fun_handler = opt_input.error_fun;
        error_fun = 1;
        error_fun_handler = error_fun_handler;
        [error_iterations, vars_prob] = error_fun_handler(X, x_opt, vars_prob);
    else
        error_iterations = norm(X{1}-x_opt)/norm(x_opt);
    end
    iter_for_errors = zeros(10,2);
    iter_for_errors(:,1) = 10.^(-(1:10))';    % The first column has the errors
    iter_for_errors(:,2) = Inf*ones(10,1);    % The second column the number of iterations
    iterations = 0;
end
% =========================================================================


% =========================================================================
% Algorithm

for k = 1 : MAX_ITER    

    X_prev_prev = X_prev;
    X_prev = X;
    for p = 1 : P
        if ~Stop(p)
            neighbs = neighbors{p};
            Dp = length(neighbs);
            
            % Determine the sum of the X's of the neighbors       
            sum_neighbs_X = zeros(n,1);
            for j = 1 : Dp
                neigh = neighbs(j);
                sum_neighbs_X = sum_neighbs_X + W(p,neigh)*X_prev{neighbs(j)};
            end
            
            Z{p} = W(p,p)*X{p} + sum_neighbs_X;
            Y{p} = par3*Z{p} + par2*X_prev{p} + par1*X_prev_prev{p};
            X{p} = weight*Y{p} + (1-weight)*Z{p};
        end
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
        iter_for_errors(intersect,2) = iterations + 2;
        
        iterations = iterations + 1;
    end
    
    % Stopping criterion for the outer loop
    Stop = Oreshkin_stopping_criterion(Stop, vars_network, X, X_prev,...
    EPS, EPS_OPT, ERROR_FUN, error_fun_handler, x_opt, vars_prob, k);
    
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



function [Stop] = Oreshkin_stopping_criterion(Stop, vars_network, X, X_prev,...
    EPS, EPS_OPT, ERROR_FUN, error_fun_handler, x_opt, vars_prob, k)


if k <= 2 
    return
end

% Network variables
P = vars_network.P;


if ~isempty(EPS_OPT)
    if ERROR_FUN == 1        
        [new_error, vars_prob] = error_fun_handler(X, x_opt, vars_prob);
    else
        new_error = norm(X{1} - x_opt)/norm(x_opt);
    end
    if new_error <= EPS_OPT
    	Stop = ones(P,1);
    end
    
    return;
end

for p = 1 : P
    if norm(X{p} - X_prev{p})/norm(X_prev{p}) <= EPS
        Stop(p) = 1;
    end
    
end

end


