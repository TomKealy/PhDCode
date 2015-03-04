function [x,vars_prob] = minimize_quad_prog_plus_l1_BB(p, v, c, X, vars_prob)

% [x,vars_prob] = minimize_quad_prog_plus_l1_BB(p, v, c, X, vars_prob)
%
% Solves the optimization problem
% 
%             minimize    (1/P)||x||_1 + v'*x + c||x||^2       (1)  
%             subject to  Ax = b
%
% where ||.||_1 is the L1-norm and ||.|| is the L2-norm. This function is
% used to solve Basis Pursuit (BP) with D-ADMM, when the matrix A_BP is
% partitioned vertically, i.e., each node has m_p = m/P rows of A_BP and 
% m_p entries of b_BP, where P is the number of nodes in the network.
%
% The inputs p, v, c, and X are explained in the D-ADMM documentation. We
% now explain the contents of the struct vars_prob.
%
% Inputs:
%     - p, v, c, X: explained in D-ADMM documentation
%     - vars_prob: is a struct with the following fields:
%            . handler: is a handler for this function (see D-ADMM doc.)
%            . A_BP: is the full BP matrix (all nodes) 
%            . b_BP: is the full BP vector b (all nodes)
%            . dual_var: is a cell of size P. Each component stores a
%                        vector of size m_p, the number of rows each node
%                        stores.
%
% Outputs:
%     - x: solution of (1)
%     - vars_prob: the same struct as the input, but with the field
%                  dual_var changed
%
%
% Method: This function solves the dual problem of (1) using a
% Barzilai-Borwein (BB) algorithm. It uses a warm-start for the dual
% variable initialization, obtained in the solution of the previous problem
% (for the same node).


% =========================================================================
% Initialization

P = length(X);                             % Number of nodes

% From struct vars_prob
A_BP = vars_prob.A_BP;
b_BP = vars_prob.b_BP;
dual_var = vars_prob.dual_var;

lambda0 = dual_var{p};
m_p = length(lambda0);
A = A_BP(1+(p-1)*m_p : p*m_p, :);
b = b_BP(1+(p-1)*m_p : p*m_p);


% Parameter of the BB algorithm
MAX_ITER = 200;
M = 20;
gamma = 1e-4;
epsilon = 1e-10;
sigma = 0.3;
alpha = 1;

n = size(A,2);
% =========================================================================


% =========================================================================
% Algorithm

lambda = lambda0;

[f_lambda, g_lambda, x] = feval(@f_and_g,lambda,A,b,v,c,n,P);


f_r = -inf*ones(M,1);
f_r(1) = f_lambda;
f_ref = f_lambda;
aux = 1;

counter = 0;

for k = 1 : MAX_ITER
    if (alpha <= epsilon) || (alpha >= 1/epsilon)
        if norm(g_lambda) > 1
            delta = 1;
        elseif (norm(g_lambda) >= 1e-5) && (norm(g_lambda) <= 1)
            delta = 1/norm(g_lambda);
        else % norm(g) < 1e-5
            delta = 1e5;
        end;
                     
        alpha = delta;
    end;
    
    eta = 1/alpha;
    
%    f_ref = max(f_r);
    if counter >= M
        f_ref = max(f_r);
        counter = 0;
    else if f_lambda > f_ref
            f_ref = f_lambda;
            counter = 0;
        end
    end
    counter = counter + 1;     
    
    slope = gamma*(g_lambda'*g_lambda);
    iter = 0;
    [f_lambda_eta_g, g_dummy, x] = feval(@f_and_g,lambda-eta*g_lambda,A,b,v,c,n,P);
    while f_lambda_eta_g > f_ref - eta*slope
        eta = sigma*eta;
        [f_lambda_eta_g, g_dummy, x] = feval(@f_and_g,lambda-eta*g_lambda,A,b,v,c,n,P);
        iter = iter + 1;
    end
    
    
    lambda_new = lambda - eta*g_lambda;
    [f_new, g_new, x] = feval(@f_and_g,lambda_new,A,b,v,c,n,P);
    
    y = g_new-g_lambda;
    
    alpha = -g_lambda'*y/(eta*(g_lambda'*g_lambda));    
    
    lambda = lambda_new;
    f_lambda = f_new;
    g_lambda = g_new;
    
    if norm(g_lambda) <= 1e-6        
        break;
    end
           
    aux = mod(aux,M)+1;
    
    f_r(aux) = f_lambda;      
    
end


end

function [f_lambda, g_lambda, x] = f_and_g(lambda,A,b,v,c,n,P)

u = v - A'*lambda;

x = zeros(n,1);

one_over_P = 1/P;

ind_pos = u < -one_over_P;
ind_neg = u > one_over_P;

x(ind_pos) = -(u(ind_pos) + one_over_P)/(2*c);
x(ind_neg) = (one_over_P - u(ind_neg))/(2*c);

g_lambda = A*x - b;

f_lambda = -(one_over_P*norm(x,1) + (u+c*x)'*x + b'*lambda);

end





