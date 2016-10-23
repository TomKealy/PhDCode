function demo_nonnegative_20110628
% 
%  YALL1_group demo of group-sparse basis pursuit 
%  with nonnegativity constraints
%
%  Recover a nonnegative group-sparse vector x by solving
%
%      minimize     sum_i w_i*||x_{gi}||_2
%      subject to   Ax = b, x>=0,
%
%  where gi is the index set of the i-th group, and
%        w_i is the weight for the i-th group.
%
%  [Adapted from spgdemo.m in the SPGL1 package.]
% -----------------------------------------------------------
clc
% Initialize random number generator
% Uncomment below to use specific seed
%     seed = 2011;
%     fprintf('seed = %d\n', seed);
%     RandStream.setDefaultStream(RandStream('mt19937ar','seed',seed));


% Set problem size
n = 1024;       % dimension of x
m = ceil(n/4);  % number of measurements

% Create random m-by-n encoding matrix
A = randn(m,n); 
% Scaling by normalizing the rows of A
d = 1./sqrt(sum(A.^2,2));
A = bsxfun(@times,A,d);
nonorth = true;  % A is non-orthonormal

% Set the numbers of groups and nonzero groups
nGroups = n / 8;
K = nGroups / 8;
groups = [];

% Generate groups with the desired number of unique groups
while (length(unique(groups)) ~= nGroups)
    groups  = ceil(rand(n,1) * nGroups);
end

% Determine weight for each group
weights = ones(nGroups,1);

% Create K-group-sparse vector x0 and observation vector b
p   = randperm(nGroups); p = p(1:K);
idx = ismember(groups,p);
nNZ = sum(idx);  % number of nonzeros
x0  = zeros(n,1); x0(idx) = abs(randn(sum(idx),1));

% Generate measurements
b   = A*x0;

% Add noise to measurements
noise_std = 0*5e-3;
noise = randn(m,1);
b     = b + noise_std*norm(b)/norm(noise)*noise;

% Set stopping tolerance value
tol = max(0.1*noise_std, 1e-6);

fprintf('--- n    m    nG    K   nNZ     std      tol ---\n')
fprintf('%5i %5i  %3i  %3i   %3i   %4.1e  %4.1e \n',...
    n,m,nGroups,K,nNZ,noise_std,tol)
fprintf('true x >= 0\n\n');

%% Call Primal YALL1_Group Solver with Exact Subproblem Solving and Nonnegative Constraints x >= 0

beta = [1,1]/mean(abs(b)); % penalty parameters
[x_yp,Out_yp] = YALL1_group(A,b,groups,...
                               'Nonnegative', true, ...     % x >= 0
                               'StopTolerance', tol, ...
                               'overlap', false, ...
                               'GrpWeights', weights, ...
                               'nonorthA', nonorth, ...
                               'ExactLinSolve', true, ...
                               'Solver', 1, ...             % 1 - primal; 2 - dual
                               'QuadPenaltyPar', beta, ...
                               'maxIter', 500);
relerr = norm(x_yp-x0)/norm(x0);  % relative error
fprintf('YALL1 Primal (Exact: true, Continuation: false, Nonnegative: true): \n')
fprintf('iter %4i, time %6.2f, relerr %6.2e\n\n',...
    Out_yp.iter,Out_yp.cputime,relerr)


%% Call Dual YALL1_Group Solver with Inexact Subproblem Solving, Continuation on the Penalty Parameters, and x >= 0

beta = 0.5*mean(abs(b)); % penalty parameter
[x_yd,Out_yd] = YALL1_group(A,b,groups,...
                               'Nonnegative', true, ...     % x >= 0
                               'StopTolerance', tol, ...
                               'GrpWeights', weights, ...
                               'overlap', false, ...
                               'nonorthA', nonorth, ...
                               'ExactLinSolve', false, ...
                               'Solver', 2, ...             % 1 - primal; 2 - dual
                               'QuadPenaltyPar', beta, ...
                               'maxIter', 500, ...
                               'Continuation', true, ...
                               'ContFactor', 1.2, ...
                               'ContParameter', 0.9);
relerr = norm(x_yd - x0) / norm(x0);
fprintf('YALL1 Dual   (Exact: false, Continuation: true, Nonnegative: true): \n')
fprintf('iter %4i, time %6.2f, relerr %6.2e\n\n',...
    Out_yd.iter,Out_yd.cputime,relerr)

%% Call SPGL1
if exist('spgl1','file')
    opts = spgSetParms('verbosity', 0, ...
                       'weights',   weights, ...
                       'bpTol',     tol, ...
                       'optTol',    5*tol, ...
                       'decTol',    tol, ...
                       'iterations',1e3);
    sigma = noise_std*norm(b);
    [x_spg,~,~,info] = spg_group(A,b,groups,sigma,opts);
    relerr = norm(x_spg - x0) / norm(x0);
    fprintf('SPGL1:\niter %4i, time %6.2f, relerr %6.2e\n',...
        info.iter,info.timeTotal,relerr)
else
    fprintf('SPGL1: not found\n')
end

end

