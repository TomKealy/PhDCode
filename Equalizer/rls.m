function [e, w] = rls(lambda, M, u, d, delta)
w = zeros(M, 1);
P = eye(M)/delta;
u = u(:) ;
d = d(:);
mu = lambda^(-1);

N = length(u);

e = d;

for n = M:N
    uvec = u(n:-1:n-M+1);
    k = mu*P*uvec/(1 + mu*uvec'*P*uvec);
    e(n) = d(n) - w'*uvec;
    w = w + k*conj(e(n));
    P = mu*P - mu*k*uvec'*P;
end
end
