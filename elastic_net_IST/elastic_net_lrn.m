function [beta,k] = elastic_net_lrn(X,Y,eps,lambda,stop,tol)
%%%%%%%%%%%%%%%%      USAGE   %%%%%%%%%%%%%%%%%
%ELASTIC_NET_LRN    iteratively solves the  elastic net regularization.
%   [BETA] = ELASTIC_NET_LRN(X,Y,EPS,LAMBDA) returns the coefficient vector
%   which solves the elastic-net regularization problem
%        min { || X BETA -Y ||^2 + LAMBDA( |BETA|^2_2 + EPS |BETA|_1} 
%   If the input data X is a NxD matrix, and the labels Y are a Nx1 vector,
%   BETA is a Dx1 vector. The solution BETA is computed via iterative 
%   soft-thresholding, with damping factor 1/(1+EPS*LAMBDA), thresholding 
%   factor EPS*LAMBDA, null initialization vector and step 
%   1/(eig_max(X*X')*1.1). The algorithm stops when the support of BETA 
%   reached convergence.
%
%   [BETA,K] = ELASTIC_NET_LRN(X,Y,EPS,LAMBDA) also returns the number of
%   iterations.
%
%   [...] = ELASTIC_NET_LRN(X,Y,EPS,LAMBDA,STOP) if STOP=0 the algorithm 
%   stops when the support of BETA reached convergence; if STOP=1 the 
%   algorithm stops when the coefficients reached convergence, that is when
%   (BETA_{l}(i)-BETA_{l+1}(i))>0.01*BETA_{l}(i) for all i.
%
%   [...] = ELASTIC_NET_LRN(X,Y,EPS,LAMBDA,STOP,TOL) if STOP=1 the 
%   algorithm stops when the coefficients reached convergence, that is when
%   the (BETA_{l}(i)-BETA_{l+1}(i))>TOL*BETA_{l}(i) for all i.
%   ELASTIC_NET_LRN(X,Y,EPS,LAMBDA) = ELASTIC_NET_LRN(X,Y,EPS,LAMBDA,0)
%   = ELASTIC_NET_LRN(X,Y,EPS,LAMBDA,0,TOL).
%   ELASTIC_NET_LRN(X,Y,EPS,LAMBDA,1) = ELASTIC_NET_LRN(X,Y,EPS,LAMBDA,1,0.01)

if nargin<4, error('too few input'), end
if nargin<5; stop = 0; end
if nargin<6; tol = 0.01; end
if nargin>6, error('too many input'), end

% start initialization
XT=X';
XY=XT*Y;
sigma_max=normest(X*XT); % largest eigenvalue
step=1/(sigma_max*1.1);
lambda = lambda*length(Y)*step;
damp_fac = 1/(1+lambda*eps); % damping factor
beta0 = zeros(length(XY),1);

k = 0; % iteration index
kmax = 100000; %the maximum number of iterations
imax = 10; % number of times the convergence condition must be satisfied
           % to stop the iteration
i = 0; % counter for the convergence condition

% first iteration:
beta = thresholding(beta0+step*(-XT*(X*beta0)+XY),lambda)*damp_fac; 

% end initialization

% iterative soft thresholding

% stopping criterion: the algorithm stops if the convergence conditions 
% below are satisfied imax times or if a maximum number of iteration kmax 
% is reached

while and(k<kmax,i<imax)
    % convergence conditions
    if stop % check convergence within a threshold (tol) for all coefficients in beta  
        if any(abs(beta0-beta) > tol*abs(beta0));
            i = 0;
        end
    else % check convergence of the support of  beta (location of non zero coefficients)
        if any(((beta0~=0)-(beta~=0))~=0);
            i = 0;
        end
    end
    beta0 = beta;
    beta = thresholding(beta0+step*(-XT*(X*beta0)+XY),lambda)*damp_fac;
    i = i+1;
    k = k+1;
    imax = max(1000,10^floor(log10(k))); %imax must be of the order of k
end

function [beta] = thresholding(beta0,lambda)
% THRESHOLDING computes the soft thresholding
%   BETA = THRESHOLDING(BETA0,LAMBDA) returns the resulting coefficient 
%   vector of the soft-thresholding of vector BETA0 with threshold
%   LAMBDA/2.

ind = logical(abs(beta0)<lambda/2);
beta = beta0-sign(beta0).*lambda/2;
beta(ind) = 0;
