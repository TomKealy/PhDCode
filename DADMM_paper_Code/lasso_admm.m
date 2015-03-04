% Function to perform LASSO regression using Alternating Direction Method
% of Multipliers.
%
% arg min_{B} 0.5*||b - A*x||_{2}^{2} + gamma*||x||_{1}
%
% Usage:- [B,cost] = lasso_admm(X, A, gamma) 
%
% where:- <in>
%         b = bias vector
%         lambda = weighting on the l1 penalty
%         <out>
%         x = solution  
%
% Written by Simon Lucey 2012

function [B,cost] = lasso_admm(b, A, gamma)

% Get dimensions of B

[m,n] = size(A)

c = size(b,2);
r = size(A,2); 

La = zeros(n,1); % Initialize Lagragian to be nothing (seems to work well)

C = zeros(n,1); % Initialize C to zero

B = zeros(n,1); 

rho = 1e-4; % Set rho to be quite low to start with 
maxIter = 500; % Set the maximum number of iterations (make really big to ensure convergence)
I = speye(r); % Set the sparse identity matrix
maxRho = 5; % Set the maximum mu


% Set the fast soft thresholding function
fast_sthresh = @(x,th) sign(x).*max(abs(x) - th,0);

% Set the norm functions
norm2 = @(x) x(:)'*x(:); 
norm1 = @(x) sum(abs(x(:))); 

cost = [];

Atb = A'*b
[L U] = factor(A, rho);

for n = 1:maxIter
        
    % Solve sub-problem to solve B
    eta = (Atb + rho*(C - La));
    nu = inv((A'*A + rho*I));
    B = eta/rho - (A'*(U \ ( L \ (A*eta) )))/(rho^2);
    
    % Solve sub-problem to solve C
    C = fast_sthresh(B + La/rho, gamma/rho); 
    
    % Update the Lagrangian
    La = La + rho*(B - C);  
    
    %pause; 
    
    % Section 3.3 in Boyd's book describes strategies for adapting rho
    % main strategy should be to ensure that 
    rho = min(maxRho, rho*1.1); 
    
    % get the current cost
    cost(n) = 0.5*norm2(b - A*B) + gamma*norm1(B);
end

end

function [L U] = factor(A, rho)
    [m, n] = size(A);
    if ( m >= n )    % if skinny
       L = chol( A'*A + rho*speye(n), 'lower' );
    else            % if fat
       L = chol( speye(m) + 1/rho*(A*A'), 'lower' );
    end
    
    % force matlab to recognize the upper / lower triangular structure
    L = sparse(L);
    U = sparse(L');
end