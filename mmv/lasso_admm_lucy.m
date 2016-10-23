% Function to perform LASSO regression using Alternating Direction Method
% of Multipliers.
%
% arg min_{B} 0.5*||X - A*B||_{2}^{2} + gamma*||B||_{1}
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

function [B,cost] = lasso_admm_lucy(X, A, gamma, rho)

% Get dimensions of B
c = size(X,2);
r = size(A,2); 

L = zeros(r,c); % Initialize Lagragian to be nothing (seems to work well)
%rho = 1e-4; % Set rho to be quite low to start with 
maxIter = 1000; % Set the maximum number of iterations (make really big to ensure convergence)
I = speye(r); % Set the sparse identity matrix
maxRho = 5; % Set the maximum mu
C = randn(r,c); % Initialize C randomly

% Set the fast soft thresholding function
fast_sthresh = @(x,th) sign(x).*max(abs(x) - th,0);

% Set the norm functions
norm2 = @(x) x(:)'*x(:); 
norm1 = @(x) sum(abs(x(:))); 

cost = [];
for n = 1:maxIter
        
    % Solve sub-problem to solve B
    B = (A'*A+rho*I)\(A'*X + rho*C - L); 
    
    % Solve sub-problem to solve C
    C = wthresh(B + (L/rho), 's', gamma/rho); 
    
    % Update the Lagrangian
    L = L + rho*(B - C);  
    
    %pause; 
    
    % Section 3.3 in Boyd's book describes strategies for adapting rho
    % main strategy should be to ensure that 
    %rho = min(maxRho, rho*1.1); 
    
    % get the current cost
    cost(n) = 0.5*norm2(X - A*B) + gamma*norm1(B);
    
end
end

function z = shrink_group(x,y,nonneg,groups,g)
%
%  Group-wise or row-wise shrinkage operator: Shrink(x,y)
%

if nonneg  % nonnegativity
    x = max(x,0);
end

t = max(0,1-y./sqrt(sum(x.^2,2)));

z = bsxfun(@times,t,x);

end