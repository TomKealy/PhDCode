% Least angle regression (LAR) algorithm
% Author: Xiaohui Chen (xiaohuic@ece.ubc.ca)
% Version: 2012-Feb

function [beta, A, mu, C, c, gamma] = lar(X, Y, option, t, standardize)

% Least Angle Regression (LAR) algorithm.
% Ref: Efron et. al. (2004) Least angle regression. Annals of Statistics.
% option = 'lar' implements the vanilla LAR algorithm (default);
% option = 'lasso' solves the lasso path with a modified LAR algorithm.
% t -- a vector of increasing positive real numbers. If given, LARS stops and 
% returns the solution at t.
%
% Output:
% A -- a sequence of indices that indicate the order of variable inclusions;
% beta: history of estimated LARS coefficients;
% mu -- history of estimated mean vector;
% C -- history of maximal current absolute corrrelations;
% c -- history of current corrrelations;
% gamma: history of LARS step size.
% Note: history is traced by rows. If t is given, beta is just the
% estimated coefficient vector at the constraint ||beta||_1 = t.
%
% Remarks:
% 1. LARS is originally proposed to estimate a sparse coefficient vector in
% a noisy over-determined linear system. LARS outputs estimates for all
% shrinkage/constraint parameters (homotopy).
%
% 2. LARS is well suited for Basis Pursuit (BP) purpose in the real case. This lars function
% automatically terminates when the current correlations for inactive set are
% all zeros. The recovered coefficient vector is the last column of beta 
% with the *lasso* option. Hence, this function provides a fast and 
% efficient solution for the ell_1 minimization problem. 
% Ref: Donoho and Tsaig (2006). Fast solution of ell_1 norm minimization problems when the solution may be sparse.

if nargin < 5, standardize = true; end
if nargin < 4, t = Inf; end
if nargin < 3, option = 'lar'; end

if strcmpi(option, 'lasso'), lasso = 1; else, lasso = 0; end

eps = 1e-10;    % Effective zero

[n,p] = size(X);
if standardize,
    X = normalize(X);
    Y = Y-mean(Y);
end
m = min(p,n-1); % Maximal number of variables in the final active set
T = length(t);

beta = zeros(1,p);
mu = zeros(n,1);    % Mean vector
gamma = []; % LARS step lengths
A = [];
Ac = 1:p;
nVars = 0;
signOK = 1;
i = 0;
mu_old = zeros(n,1);
t_prev = 0;
beta_t = zeros(T,p);
ii = 1;
tt = t;

% LARS loop
while nVars < m,
    i = i+1;
    c = X'*(Y-mu);  % Current correlation
    C = max(abs(c));    % Maximal current absolute correlation
    if C < eps || isempty(t), break; end    % Early stopping criteria
    if 1 == i, addVar = find(C==abs(c)); end
    if signOK,
        A = [A,addVar]; % Add one variable to active set
        nVars = nVars+1;
    end
    s_A = sign(c(A));
    Ac = setdiff(1:p,A);    % Inactive set
    nZeros = length(Ac);
    X_A = X(:,A);
    G_A = X_A'*X_A; % Gram matrix
    invG_A = inv(G_A);
    L_A = 1/sqrt(s_A'*invG_A*s_A);
    w_A = L_A*invG_A*s_A;   % Coefficients of equiangular vector u_A
    u_A = X_A*w_A;  % Equiangular vector
    a = X'*u_A; % Angles between x_j and u_A
    beta_tmp = zeros(p,1);
    gammaTest = zeros(nZeros,2);
    if nVars == m,
        gamma(i) = C/L_A;   % Move to the least squares projection
    else
        for j = 1:nZeros,
            jj = Ac(j);
            gammaTest(j,:) = [(C-c(jj))/(L_A-a(jj)), (C+c(jj))/(L_A+a(jj))];
        end
        [gamma(i) min_i min_j] = minplus(gammaTest);
        addVar = unique(Ac(min_i));
    end
    beta_tmp(A) = beta(i,A)' + gamma(i)*w_A;    % Update coefficient estimates
    % Check the sign feasibility of lasso
    if lasso,
        signOK = 1;
        gammaTest = -beta(i,A)'./w_A;
        [gamma2 min_i min_j] = minplus(gammaTest);
        if gamma2 < gamma(i),   % The case when sign consistency gets violated
            gamma(i) = gamma2;
            beta_tmp(A) = beta(i,A)' + gamma(i)*w_A;    % Correct the coefficients
            beta_tmp(A(unique(min_i))) = 0;
            A(unique(min_i)) = [];  % Delete the zero-crossing variable (keep the ordering)
            nVars = nVars-1;
            signOK = 0;
        end
    end
    if Inf ~= t(1),
        t_now = norm(beta_tmp(A),1);
        if t_prev < t(1) && t_now >= t(1),
            beta_t(ii,A) = beta(i,A) + L_A*(t(1)-t_prev)*w_A';    % Compute coefficient estimates corresponding to a specific t
            t(1) = [];
            ii = ii+1;
        end
        t_prev = t_now;
    end
    mu = mu_old + gamma(i)*u_A; % Update mean vector
    mu_old = mu;
    beta = [beta; beta_tmp'];
end

if 1 < ii,
    noCons = (tt > norm(beta_tmp,1));
    if 0 < sum(noCons),
        beta_t(noCons,:) = repmat(beta_tmp',sum(noCons),1);
    end
    beta = beta_t;
end


% Normalize columns of X to have mean zero and length one.
function sX = normalize(X)

[n,p] = size(X);
sX = X-repmat(mean(X),n,1);
sX = sX*diag(1./sqrt(ones(1,n)*sX.^2));


% Find the minimum and its index over the (strictly) positive part of X
% matrix
function [m, I, J] = minplus(X)

% Remove complex elements and reset to Inf
[I,J] = find(0~=imag(X));
for i = 1:length(I),
    X(I(i),J(i)) = Inf;
end

X(X<=0) = Inf;
m = min(min(X));
[I,J] = find(X==m);
