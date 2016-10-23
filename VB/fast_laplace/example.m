% This is a simple example to test the fast Laplace algorithm from the following paper:
%
% [1] S. D. Babacan, R. Molina, A. K. Katsaggelos. “Bayesian Compressive Sensing using Laplace Priors,” 
% submitted for publication, IEEE Transactions on Image Processing, September 2008.
% 
% Author: S. Derin Babacan

M = 100; % number of measurements
N = 512; % length of the signal
T = 20; % number of nonzero components

%% Generate the sparse signal
randn('state',sum(100*clock))
rand('state',sum(100*clock))

% Random +-1 spikes
q = randperm(N);
x = zeros(N,1);
x(q(1:T)) = sign(randn(T,1));

A = rand(M,N);
A = A./repmat(sqrt((sum(A.^2,1))),[M,1]);

%% Generate the observations

% Observation noise standard deviation
sigma = 0.005; 

e = sigma*randn(M,1);
y = A*x + e;

%% Run the reconstruction
initsigma2 = std(y)^2/1e2;

% lambda = 0; % Gaussian priors (BCS)
% lambda = <a_positive_value> -- Laplace priors with fixed lambda
lambda = []; % Laplace priors with lambda calculated adaptively (see [1])

tic;
[weights,used,sigma2,errbars,basis,selected,alpha,lambdas] = FastLaplace(A,y,initsigma2,1e-8,lambda);
t_Lap = toc;
x_Lap = zeros(N,1);  x_Lap(used) = weights;
err = zeros(N,1); err(used) = errbars;
s_Lap = length(used);
E_Lap = norm(x-x_Lap)/norm(x);

%% Show the results
fprintf(1,'Number of nonzero weights: %d,  time = %g, err = %g\n',length(used), t_Lap, E_Lap);

% Show the reconstruction
figure(1), subplot(2,1,1), stem(x), title('Original Signal'), axis([1 N -2 2]);
subplot(2,1,2),stem(x_Lap), axis([1 N -2 2]);, title(sprintf('Reconstruction from %d measurements, error = %g',M,E_Lap));

% Show the reconstruction with estimation uncertainties
figure(2), subplot(2,1,1), stem(x), title('Original Signal'), axis([1 N -2 2]);
subplot(2,1,2),errorbar(x_Lap,err);, axis([1 N -2 2]);, title(sprintf('Reconstruction from %d measurements with error bars',M));





            