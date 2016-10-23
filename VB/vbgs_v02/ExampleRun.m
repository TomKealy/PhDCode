% Example run of VBGS.m
% Variational Bayesian Group-Sparse Estimation
%    S. D. Babacan, M. Luessi, R. Molina, and A. K. Katsaggelos, 
%    "Sparse Bayesian Methods for Low-Rank Matrix Estimation," 
%    IEEE Transactions on Signal Processing, 2012.
%
% 
clear
clc
close all

N = 500; % length of vector w
NoNonZero= 100; % no of nonzeros in the unknown vector w
group_size = 20; 
NoNonZeroGroups = NoNonZero / group_size;
nogroups = N/group_size;
M = 160; % no measurements


%% Generate random dictionary (measurement matrix)
A = randn(M,N);
A = A./repmat(sqrt(sum(A.^2,1)),[M,1]);

%% Generate signal
w = zeros(N,1);
%Grouping: Simple sequential block structure, but can be arbitrary too
groups = reshape( ([1:nogroups]'*ones(1,group_size))' ,[N,1]);

for i=1:nogroups,
    grouping{i} = find(ismember(groups, i));
end

% Generate nonzero groups
ind = randperm(nogroups); ind = ind(1:NoNonZeroGroups);
ind = ismember(groups, ind);
w(ind) = randn(NoNonZero,1);

% Generate observation
y = A*w;
sigma = 1e-3;% some noise
% sigma = sqrt(1e-3); % more noise
y = y + sigma*randn(M,1); 

%% OPTIONS
% All options are optional, i.e., you can put options=[] and default 
% options will be generated.
% see VBGS.m for explanations
% These need not be modified unless testing method behaviour
options.verbose              = 0; % verbosity
options.lambda               = 1; 
options.k_a                  = 1e-6; % Hyperparameters
options.theta_a              = 1e-6; 
options.k_b                  = 1e-6;
options.theta_b              = 1e-6;
options.k_beta0              = 1e-6;
options.theta_beta0          = 1e-6;
options.MAXITER              = 200; 
options.pt_estimate          = 'Mean';
options.conv_thr             = 1e-7; % convergence threshold
options.ESTIMATE_HYPERPARAMS = 1; % estimate hyperparameters?
% If the hyperparameters a,b are NOT to be estimated, supply a0 and b0 below
% This is not recommended, updating them via Bayesian estimation generally
% gives much better results
options.a0                   = 1;
options.b0                   = 1; 


%% The following can be modified 
% estimate noise variance? 
options.UPDATE_BETA          = 1; 
options.beta_init            = 'auto'; % Automatic
%Otherwise set it to a value (for instance, 1/sigma), 
%and uncomment 2 lines below.
% options.UPDATE_BETA = 0
% options.beta_init            = 1/sigma;
% tuning the noise precision can give superior results compared to 
% automatically estimating it. use small beta_init for high noise, high
% beta_init for low noise.

options.UPDATE_BETA_init     = 1; % iteration number to start estimating the noise variance,
% if the magnitude difference within the signal is too large, set it higher
% (like 3)

% Covariance estimation
% Set to 2 for small to medium scale, 3 for large scale problems.
% Set 4 for MAP solution (covariances 0), 1 is for standard VB (slow, see
% paper)
options.cov_met              = 3; 


% Dimensionality reduction (not included in the paper). 
% This is basically a thresholding operation on the calculated inverse 
% variances at each iteration. It can significantly decrease the running
% time. It currently does not work with overlapping groups.
options.DIMRED  = 1;
% You can also set the threshold via
options.DIMRED_THR = 1e2;

% signal distribution
% Choose from 'McKay', 'MLaplace', 'MStudent' or 'Jeffreys'
% Default : 'Jeffreys'
% options.dist                 = 'McKay'; 
% options.dist                 = 'MStudent'; 
% options.dist                 = 'MLaplace'; 
options.dist                   = 'Jeffreys'; 




%% Solve using VBGS
tic 
[w_est, it, z, beta, Sigma_w, a, b, used] = VBGS(y, A, grouping, options);
t = toc;
err = norm( w_est - w , 'fro') / norm( w , 'fro');

%% Show results
fprintf(1,'VBGS example run \n')
fprintf(1,'Signal length = %d, group size = %d, %d nonzero groups\n', N, group_size, NoNonZeroGroups);
fprintf(1,'%d observations, noise variance %g\n', M, sigma)
if ~isempty(options) & isfield(options, 'dist'),
    fprintf(1,'Signal prior = %s\n', options.dist); 
else
    fprintf(1,'Signal prior = Jeffreys\n'); 
end
fprintf(1,'Relative Estimation error = %g, computation time = %g s\n', err, t);

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])

subplot(1,3,1), stem([1:M], y)
axis([1, M, min(y(:)), max(y(:))])
th = title('Observation');
set([gca th],'FontSize', 20)

subplot(1,3,2),stem([1:N], w), 
axis([1, N, min(w(:)), max(w(:))])
th = title('Original');
set([gca th],'FontSize', 20)

subplot(1,3,3), stem([1:N], w_est), 
axis([1, N, min(w_est(:)), max(w_est(:))])
th = title('Estimated');
xlh = xlabel(sprintf('Rec. err. =%.3f', err));
set([gca th xlh],'FontSize', 20)
