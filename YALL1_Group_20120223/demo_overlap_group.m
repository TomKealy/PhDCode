function demo_overlap_group_20120219
% 
%  YALL1_group demo of overlapping group-sparse basis pursuit
%
%  Recover a group-sparse vector x by solving
%
%      minimize     sum_i w_i*||x_{gi}||_2
%      subject to   Ax = b,
%
%  where gi is the index set of the i-th group, w_i is the weight
%  for the i-th group, and the groups may overlap.
%
% -------------------------------------------------------------------
clc
% Initialize random number generator
% Uncomment below to use specific seed
%     seed = 2011;
%     fprintf('seed = %d\n', seed);
%     RandStream.setDefaultStream(RandStream('mt19937ar','seed',seed));


% Set problem size
n = 1024;       % dimension of x
m = ceil(n/3);  % number of measurements

% Create random m-by-n encoding matrix
A = randn(m,n); 
% Scaling by normalizing the rows of A
d = 1./sqrt(sum(A.^2,2));
A = bsxfun(@times,A,d);

% Set the numbers of groups and nonzero groups
nGroups = n / 8;
K = nGroups / 8;
groups = cell(nGroups,1);

% Generate overlapping groups
for i = 1:nGroups-1
    groups{i}  = 8*(i-1)+(1:10)';
end
groups{nGroups} = (n-7:n)';

% Determine weight for each group
weights = ones(nGroups,1);

% Create K-group-sparse vector x0 and observation vector b
p   = randperm(nGroups);
x0  = randn(n,1); 
for i = 1:nGroups-K
    x0(groups{p(i)}) = 0;
end

% Generate measurements
b   = A*x0;

% Add noise to measurements
noise_std = 0*5e-3;
noise = randn(m,1);
b     = b + noise_std*norm(b)/norm(noise)*noise;

% Set stopping tolerance value
tol = max(0.1*noise_std, 1e-6);

fprintf('--- n    m    nG    K      std      tol ---\n')
fprintf('%5i %5i   %3i  %3i   %4.1e  %4.1e \n\n',...
    n,m,nGroups,K,noise_std,tol)

%% Call Primal YALL1_Group Solver with Inexact Subproblem Solving

beta = 0.1*[1,1]/mean(abs(b)); % penalty parameters
[x_yp,Out_yp] = YALL1_group(A,b,groups,...
                               'StopTolerance', tol, ...
                               'GrpWeights', weights, ...
                               'overlap', true, ...
                               'Solver', 1, ...             % 1 - primal; 2 - dual
                               'QuadPenaltyPar', beta, ...
                               'maxIter', 500);
relerr = norm(x_yp-x0)/norm(x0);  % relative error
fprintf('YALL1 Primal (Exact: false, Continuation: false): \n')
fprintf('iter %4i, time %6.2f, relerr %6.2e\n\n',...
    Out_yp.iter,Out_yp.cputime,relerr)

%% Call Primal YALL1_Group Solver with Inexact Subproblem Solving and Continuation on the Penalty Parameters

beta = 0.01*[1,1]/mean(abs(b)); % penalty parameters
[x_yp,Out_yp] = YALL1_group(A,b,groups,...
                               'StopTolerance', tol, ...
                               'GrpWeights', weights, ...
                               'overlap', true, ...
                               'Solver', 1, ...             % 1 - primal; 2 - dual
                               'QuadPenaltyPar', beta, ...
                               'maxIter', 500, ...
                               'Continuation', true, ...
                               'ContFactor', 1.1, ...
                               'ContParameter', 0.9);
relerr = norm(x_yp - x0) / norm(x0);
fprintf('YALL1 Primal (Exact: false, Continuation: true): \n')
fprintf('iter %4i, time %6.2f, relerr %6.2e\n\n',...
    Out_yp.iter,Out_yp.cputime,relerr)


end

