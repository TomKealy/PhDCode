function alpha = SolveFastCS(alpha0, N)
% SolveFastCS: Solves a Compressed Sensing problem using fast operators.
%  Usage:
%    alpha = SolveFastCS(alpha0, N)
%  Inputs:
%    alpha0  Original data vector, of length M
%    N       number of sensed samples to use
%  Outputs:
%    alpha   recontructed vector 
%  Description:
%    SolveFastCS simulates a fast compressed sensing scheme as follows:
%    It generates a random sampling matrix Phi, size NxM, with M =
%    length(alpha0), and N specified. 
%    It then solves the BP problem
%      min || alpha ||_1 s.t. Phi*alpha = Phi*alpha0
%    The matrix Phi is constructed as P*H*Q, where P,Q are random permutation
%    matrices, and H is a fast Hadamard/Fourier operator.
%
%    The program relies on Michael Saunders' PDCO routine 
%    (primal-dual log barrier interior point method) to 
%    solve this problem, available at:
%    http://www.stanford.edu/group/SOL/software/pdco.html

% Global variables to FastCSMatrix:
% Pstate - state of random generator for P matrix
% Qstate - state of random generator for Q matrix
global Pstate Qstate

Pstate = 4972169;
Qstate = 7256157;

alpha0 = alpha0(:);
M = length(alpha0);

n = 2*M;    % Input size
m = N;      % Output size

% generate the vector b
b = FastCSMatrix(1,m,n,[alpha0; zeros(size(alpha0))]);

% generate upper and lower bounds
bl = zeros(n,1);
bu = Inf .* ones(n,1);

% generate the vector c
c = ones(n,1);

% Generate an initial guess
x0    = ones(n,1)/n;       % Initial x
y0    = zeros(m,1);        % Initial y
z0    = ones(n,1)/n;       % Initial z

d1 = 1e-4;                 % Regularization parameters
d2 = 1e-4;
xsize = 1;                 % Estimate of norm(x,inf) at solution
zsize = 1;                 % Estimate of norm(z,inf) at solution

CGMaxIter = n * log2(n) / m;
options = pdcoSet;         % Option set for the function pdco
options = pdcoSet( options, ...
                     'MaxIter    ', 50        , ...
                     'FeaTol     ', 1e-2      , ...
                     'OptTol     ', 1e-2      , ...
                     'StepTol    ', 0.99      , ...
                     'StepSame   ', 0         , ...
                     'x0min      ', 0.1       , ...
                     'z0min      ', 1.0       , ...
                     'mu0        ', 0.01      , ...
                     'method     ', 1         , ...
                     'LSQRMaxIter', 20 , ...
                     'LSQRatol1  ', 1e-3      , ...
                     'LSQRatol2  ', 1e-15     , ... 
                     'wait       ', 0    );

[x,y,z,inform,PDitns,CGitns,time] = ...
    pdco(c, 'FastCSMatrix', b, bl, bu, d1, d2, options, x0, y0, z0, xsize, zsize);

% Extract the solution from the output vector x
alpha = x(1:M) - x((M+1):(2*M));
