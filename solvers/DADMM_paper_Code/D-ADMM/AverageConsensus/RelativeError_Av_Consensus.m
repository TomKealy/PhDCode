function [error, vars_prob] = RelativeError_Av_Consensus(X, x_opt, vars_prob)

P = length(X);

x_estimate = zeros(P,1);

for p = 1 : P
    x_estimate(p) = X{p};
end

error = norm(x_estimate - x_opt)/norm(x_opt);


