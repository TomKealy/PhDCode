function [X, vars_prob, varargout] = DADMM_lasso(n, vars_prob, vars_network,...
    varargin)

OPT_OUTPUT = 0;    % 0 if there is no optional output; 1 otherwise

% =========================================================================
% Check for input and output errors
optargs = nargin - 3;

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
nout = max(nargout,1)-2;

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
T = cell(P,1);
                                 % Cell (one/node) with the dual variable
Diff = cell(P,1);                % Cell (one/node) with
                                 %    Dp*X{p} - sum_{j in neighbs}X{j}
X = cell(P,1);                   % Cell (one/node) with current estimates
X_prev = cell(P,1);              % Cell (one/node) with previous estimates

Z = cell(P,1);
Z_prev = cell(P,1);

Ap = cell(P,1);
bp = cell(P,1);
Atb = cell(P,1);

for p = 1 : P
    X{p} = zeros(n,1);
    X_prev{p} = X{p};
    Z{p} = zeros(n,1);
    Z_prev{p} = Z{p};
    U{p} = zeros(n,1);
    Diff{p} = zeros(n,1);
    T{p} = zeros(n,1);
    Ap{p} = vars_prob.A_BPDN(p,:);
    bp{p} = vars_prob.b_BPDN(p,:);
    Atb{p} = Ap{p}'*bp{p};
end

Stop = zeros(P,1);               % if Stop(p) == 1, node p has terminated


if EXISTS_X_OPT == 1
    if ERROR_FUN == 1
        error_fun_handler = opt_input.error_fun;
        vars.error_fun = 1;
        vars.error_fun_handler = error_fun_handler;
        [error_iterations, vars_prob] = error_fun_handler(Z, x_opt, vars_prob);
    else
        size(x_opt)
        size(Z{1})
        error_iterations = norm(Z{1}-x_opt)/norm(x_opt);
    end
    iter_for_errors = zeros(10,2);
    iter_for_errors(:,1) = 10.^(-(1:10))';    % The first column has the errors
    iter_for_errors(:,2) = Inf*ones(10,1);    % The second column the number of iterations
    iterations = 0;
    
    vars.x_opt = x_opt;
    vars.error_iterations = error_iterations;
    vars.iter_for_errors = iter_for_errors;
    vars.iterations = iterations;
end
% =========================================================================

% =========================================================================
% Algorithm

for k = 1 : MAX_ITER
    k
    vars.X_prev = X;
    vars.Z_prev = Z;
    
    % Network variables
    P = vars_network.P;
    neighbors = vars_network.neighbors;
    partition_colors = vars_network.partition_colors;
    
    n = length(X{1});
    
    num_colors = length(partition_colors);
    
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
                    sum_neighbs = sum_neighbs + Z{neighbs(j)};
                end
                
                % Solving each problem in the alternating direction minimization
                v = U{p} - rho*sum_neighbs;
                t = T{p};
                [r,Nc] = size(partition_colors);
                zp = Z_aux{p};
                
                %x-update
                q = (Ap{p}'*Ap{p}+rho*eye(200));
                
                xp = (Atb{p}+rho*zp - t);
                
                xp = q\xp;
                
                %z-update
                zp = shrinkage(xp+((t-v)/(rho)), (vars_prob.beta/P)/(rho));
                
                X_aux{p} = xp;
                Z_aux{p} = zp;
            end
        end
        X = X_aux;
        Z = Z_aux;
    
    end
        
    for p = 1 : P
        neighbs = neighbors{p};
        Dp = length(neighbs);
        
        % Determine the sum of the X's of the neighbors
        sum_neighbs = zeros(n,1);
        for j = 1 : Dp
            sum_neighbs = sum_neighbs + Z{neighbs(j)};
        end
        
        Diff{p} = Dp*(Z{p}) - sum_neighbs;
    end
    
    % =========================================================================
    % Update Multipliers
    
    for p = 1 : P
        % Only iterate if the nodes are still active
        if Stop(p) == 0
            U{p} = U{p} + rho*Diff{p} + rho*(X{p}-Z{p});
            T{p} = T{p} + rho*(X{p}-Z{p});
        end
    end
    
    if EXISTS_X_OPT == 1
              
        if ~isempty(EPS_OPT)
               new_error = norm(Z{1} - x_opt)/norm(x_opt);
            if new_error <= EPS_OPT
                Stop = ones(P,1);
            end
            
            return;
        end
    end

for p = 1 : P
    if norm(Z{p} - Z_prev{p})/norm(Z_prev{p}) <= EPS
        Stop(p) = 1;
    end
    
end

end

function [Stop] = DADMM_stopping_criterion(Stop, vars_network, vars, vars_prob, k)


if k <= 2 || vars.TURN_OFF_EPS == 1
    return
end

% Network variables
P = vars_network.P;

% D-ADMM variables
X = vars.X;
Z = vars.Z;
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
    if norm(Z{p} - Z_prev{p})/norm(Z_prev{p}) <= EPS
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
