function y = FastCSMatrix(mode,m,n,x)
% FastCSMatrix: The operator form of the contraint matrix for the 
% fast compressed sensing problem.
% Specifically, it returns y = A*x (mode = 1) or y = A'*x (mode = 2),
% where A is an mxn matrix given by A = [Phi - Phi], and
% Phi = P*H*Q, where P,Q are random permutation matrices, 
% and H is a fast Hadamard/Fourier operator.

% Global variables:
% Pstate - state of random generator for P matrix
% Qstate - state of random generator for Q matrix
global Pstate Qstate

n2 = n/2;

if (mode == 1) % Direct operator
    % Decompose input
    u = x(1:n2);
    v = x(n2+1:n);
    
    % Apply matrix Q
    rand('state', Qstate);
    q = randperm(n2);
    u2 = u(q);
    v2 = v(q);
    
    % Apply matrix H
    u3 = RST(u2);
    v3 = RST(v2);
    
    % Apply matrix P
    rand('state', Pstate);
    p = randperm(n2);
    y = u3(p(1:m)) - v3(p(1:m));
    
else % Adjoint operator
    % Apply matrix P^T
    rand('state', Pstate);
    p = randperm(n2);
    x2 = zeros(n2,1);
    x2(p(1:m)) = x;

    % Apply matrix H^T
    x3 = Inv_RST(x2);
    
    % Apply matrix Q^T
    rand('state', Qstate);
    q = randperm(n2);
    y = zeros(n2,1);
    y(q) = x3;
    y = [y; -y];
end

