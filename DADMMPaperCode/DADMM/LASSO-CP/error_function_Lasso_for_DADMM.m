function [error, vars_prob] = error_function_Lasso_for_DADMM(X, x_opt, vars_prob)

P = length(X);

est = vars_prob.x_opt;
n_p = vars_prob.n;

n = P*n_p;

x_estimate = zeros(n,1);

for p = 1 : P
    
    x_estimate( (p-1)*n_p + 1 : p*n_p ) = est{p};
    
end

error = norm(x_estimate - x_opt)/norm(x_opt);

% % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% % Erase
% error_along_iter = vars_prob.error_along_iter;
% error_along_iter = [error_along_iter ; error];
% vars_prob.error_along_iter = error_along_iter;
% 
% figure(1);clf;
% semilogy(error_along_iter);
% drawnow;
% % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


