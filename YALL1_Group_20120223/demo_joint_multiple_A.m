function demo_joint_multiple_A
% 
%  YALL1_group demo of jointly-sparse basis pursuit
%
%  Recover a matrix X with a few nonzero rows by solving
%
%      minimize     sum_i w_i*||X(i,:)||_2
%      subject to   A{l}X(:,l) = B(:,l), for l = 1,2,...,L
%
%  where X(i,:) is the i-th row of matrix X,
%        w_i is the weight for the i-th row, and
%        A{l}X(:,l) = B(:,l) is the sensing equation of channel l.
%
% ----------------------------------------------------------------------
clc
% Initialize random number generator
% Uncomment below to use specific seed
%     seed = 2011;
%     fprintf('seed = %d\n', seed);
%     RandStream.setDefaultStream(RandStream('mt19937ar','seed',seed));

% Set problem size
n = 1024;       % dimension of x
m = ceil(n/4);  % number of measurements
L = ceil(m/16); % number of columns of X
k = 110;        % number of nonzero rows

% Create L random m-by-n encoding matrices
A = cell(L,1);
for l = 1:L
    A{l} = randn(m,n); 
    % Scaling by normalizing the rows of A
    d = 1./sqrt(sum(A{l}.^2,2));
    A{l} = bsxfun(@times,A{l},d);
end
nonorth = true;  % A is non-orthonormal

% Determine weight for each row
weights = ones(n,1);

% Create k-jointly-sparse matrix X0 and observation vector B
p = randperm(n); p = p(1:k);
X0= zeros(n,L); X0(p,:) = randn(k,L);

% Generate measurements
B = zeros(m,L);
for l = 1:L % parfor ready
    B(:,l) = A{l}*X0(:,l);
end

% Add noise to measurements
noise_std = 0*5e-3;
Noise = randn(m,L);
B = B + noise_std*norm(B(:))/norm(Noise(:))*Noise;

% Set stopping tolerance value
tol = max(0.1*noise_std, 1e-6);

fprintf('--- n    m    k    L     std      tol ---\n')
fprintf('%5i %5i  %3i  %3i  %4.1e  %4.1e \n\n',...
    n,m,k,L,noise_std,tol)


%% Call Primal YALL1_Group Solver with Exact Subproblem Solving

beta = [1,1]/mean(abs(B(:))); % penalty parameters
[X_yp,Out_yp] = YALL1_group(A,B,[],...
                               'StopTolerance', tol, ...
                               'GrpWeights', weights, ...
                               'nonorthA', nonorth, ...
                               'ExactLinSolve', true, ...
                               'Solver', 1, ...             % 1 - primal; 2 - dual
                               'QuadPenaltyPar', beta, ...
                               'maxIter', 1000);
relerr = norm(X_yp-X0,'fro')/norm(X0,'fro');  % relative error
fprintf('YALL1 Primal (Exact: true, Continuation: false): \n')
fprintf('iter %4i, time %6.2f, relerr %6.2e\n\n',...
    Out_yp.iter,Out_yp.cputime,relerr)

%% Call Primal YALL1_Group Solver with Inexact Subproblem Solving

beta = [1,1]/mean(abs(B(:))); % penalty parameters
[X_yp,Out_yp] = YALL1_group(A,B,[],...
                               'StopTolerance', tol, ...
                               'GrpWeights', weights, ...
                               'nonorthA', nonorth, ...
                               'ExactLinSolve', false, ...
                               'Solver', 1, ...             % 1 - primal; 2 - dual
                               'QuadPenaltyPar', beta, ...
                               'maxIter', 1000);
relerr = norm(X_yp-X0,'fro')/norm(X0,'fro');  % relative error
fprintf('YALL1 Primal (Exact: false, Continuation: false): \n')
fprintf('iter %4i, time %6.2f, relerr %6.2e\n\n',...
    Out_yp.iter,Out_yp.cputime,relerr)

%% Call Primal YALL1_Group Solver with Inexact Subproblem Solving and Continuation on the Penalty Parameters

beta = 0.05*[1,1]/mean(abs(B(:))); % penalty parameters
[X_yp,Out_yp] = YALL1_group(A,B,[],...
                               'StopTolerance', tol, ...
                               'GrpWeights', weights, ...
                               'nonorthA', nonorth, ...
                               'ExactLinSolve', false, ...
                               'Solver', 1, ...             % 1 - primal; 2 - dual
                               'QuadPenaltyPar', beta, ...
                               'maxIter', 1000, ...
                               'Continuation', true, ...
                               'ContFactor', 1.1, ...
                               'ContParameter', 0.9);
relerr = norm(X_yp-X0,'fro')/norm(X0,'fro');
fprintf('YALL1 Primal (Exact: false, Continuation: true): \n')
fprintf('iter %4i, time %6.2f, relerr %6.2e\n\n',...
    Out_yp.iter,Out_yp.cputime,relerr)



%% Call Dual YALL1_Group Solver with Exact Subproblem Solving

beta=1*mean(abs(B(:)));
[X_yd,Out_yd] = YALL1_group(A,B,[],...
                               'StopTolerance', tol, ...
                               'GrpWeights', weights, ...
                               'nonorthA', nonorth, ...
                               'ExactLinSolve', true, ...
                               'Solver', 2, ...             % 1 - primal; 2 - dual
                               'QuadPenaltyPar', beta, ...
                               'maxIter', 1000);
relerr = norm(X_yd-X0,'fro')/norm(X0,'fro');
fprintf('YALL1 Dual   (Exact: true, Continuation: false): \n')
fprintf('iter %4i, time %6.2f, relerr %6.2e\n\n',...
    Out_yd.iter,Out_yd.cputime,relerr)


%% Call Dual YALL1_Group Solver with Inexact Subproblem Solving

beta = 1*mean(abs(B(:))); % penalty parameters
[X_yd,Out_yd] = YALL1_group(A,B,[],...
                               'StopTolerance', tol, ...
                               'GrpWeights', weights, ...
                               'nonorthA', nonorth, ...
                               'ExactLinSolve', false, ...
                               'Solver', 2, ...             % 1 - primal; 2 - dual
                               'QuadPenaltyPar', beta, ...
                               'maxIter', 1000);
relerr = norm(X_yd-X0,'fro')/norm(X0,'fro');  % relative error
fprintf('YALL1 Dual   (Exact: false, Continuation: false): \n')
fprintf('iter %4i, time %6.2f, relerr %6.2e\n\n',...
    Out_yd.iter,Out_yd.cputime,relerr)

%% Call Dual YALL1_Group Solver with Inexact Subproblem Solving and Continuation on the Penalty Parameters

beta = 0.5*mean(abs(B(:))); % penalty parameter
[X_yd,Out_yd] = YALL1_group(A,B,[],...
                               'StopTolerance', tol, ...
                               'GrpWeights', weights, ...
                               'nonorthA', nonorth, ...
                               'ExactLinSolve', false, ...
                               'Solver', 2, ...             % 1 - primal; 2 - dual
                               'QuadPenaltyPar', beta, ...
                               'maxIter', 1000, ...
                               'Continuation', true, ...
                               'ContFactor', 1.2, ...
                               'ContParameter', 0.9);
relerr = norm(X_yd-X0,'fro')/norm(X0,'fro');  % relative error
fprintf('YALL1 Dual   (Exact: false, Continuation: true): \n')
fprintf('iter %4i, time %6.2f, relerr %6.2e\n\n',...
    Out_yd.iter,Out_yd.cputime,relerr)

end
