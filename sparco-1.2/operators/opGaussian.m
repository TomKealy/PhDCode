function op = opGaussian(m,n,mode)
% OPGAUSSIAN Gaussian ensemble
%
%    OPGAUSSIAN(M,N,MODE) creates an M by N Gaussian ensemble
%    operator. By setting MODE a number of different types of
%    ensemble can be generated; MODE = 0 (default) explicitly
%    creates an unnormalized matrix from the Normal
%    distribution. MODE = 1 generates columns of the unnormalized
%    matrix as the operator is applied. This allows for much larger
%    ensembles since the matrix is implicit. For MODE = 2,3 columns
%    are scaled to have unit Euclidean norm, when MODE = 3 the
%    matrix is implicit. For MODE 4,5 the rows are orthogonalized
%    using respectively the QR, and ORTH functions. For the last
%    two modes matrices are always generated explicitly. The
%    orthogonalization modes work only when M < N.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: opGaussian.m 1052 2008-07-02 20:52:08Z ewout78 $

if nargin < 3, mode = 0; end;


switch (mode)
 case 0
    A  = randn(m,n);
    op = opMatrix(A,{'Gaussian'});

 case 1
    seed = randn('state');
    for i=1:m, randn(n,1); end; % Ensure random state is advanced
    op   = @(x,mode) opGaussian_intrnl(m,n,seed,x,mode);

 case 2
    A  = randn(m,n);
    A  = A * spdiags((1./sqrt(sum(A.*A)))',0,n,n);
    op = opMatrix(A,{'Gaussian'});

 case 3
    seed = randn('state');
    scale = zeros(1,n);
    for i=1:n
       v = randn(m,1);
       scale(i) = 1 / sqrt(v'*v);
    end
    op = @(x,mode) opGaussianScaled_intrnl(m,n,seed,scale,x,mode);

 case 4
    if (m > n) error('This mode is not supported when M > N'); end;
    A  = randn(n,m);
    [Q,R] = qr(A,0);
    op = opMatrix(Q',{'Gaussian'});

 case 5
    if (m > n) error('This mode is not supported when M > N'); end;
    A  = randn(m,n);
    A  = orth(A')';
    op = opMatrix(A,{'Gaussian'});
    
 otherwise
     error('Invalid mode')
end



function y = opGaussian_intrnl(m,n,seed,x,mode)
checkDimensions(m,n,x,mode);

% Store current random number generator state
seed0 = randn('state');
randn('state',seed);

if mode == 0
   y = {m,n,[0,1,0,1],{'Gaussian','Unnormalized'}};
elseif mode == 1
   y = zeros(m,1);
   for i=1:n
      y = y + randn(m,1) * x(i);
   end
else
   y = zeros(n,1);
   for i=1:n
      y(i) = randn(1,m) * x;
   end
end

% Restore original random number generator state
randn('state',seed0);



function y = opGaussianScaled_intrnl(m,n,seed,scale,x,mode)
checkDimensions(m,n,x,mode);

% Store current random number generator state
seed0 = randn('state');
randn('state',seed);

if mode == 0
   y = {m,n,[0,1,0,1],{'Gaussian','Normalized'}};
elseif mode == 1
   y = zeros(m,1);
   for i=1:n
      y = y + randn(m,1) * (scale(i) * x(i));
   end
else
   y = zeros(n,1);
   for i=1:n
      y(i) = scale(i) * randn(1,m) * x;
   end
end

% Restore original random number generator state
randn('state',seed0);
