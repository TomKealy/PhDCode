function op = opSign(m,n,mode)
% OPSIGN   Sign-ensemble operator.
%
%    OPSIGN(M,N,MODE) creates an M by N sign-ensemble operator. By
%    setting MODE a number of different types of ensemble can be
%    generated; MODE = 0 (default) explicitly creates an
%    unnormalized matrix with random +1 and -1 entries. MODE = 1
%    generates columns of the unnormalized matrix as the operator
%    is applied. This allows for much larger ensembles since the
%    matrix is implicit. For MODE = 2,3 columns are scaled to have
%    unit Euclidean norm, when MODE = 3 the matrix is implicit.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: opSign.m 1040 2008-06-26 20:29:02Z ewout78 $

if nargin < 3, mode = 0; end;


switch (mode)
 case 0
    A  = 2.0 * (randn(m,n) < 0) - 1;
    op = opMatrix(A,{'Sign'});

 case 1
    seed = randn('state');
    for i=1:m, randn(n,1); end; % Ensure random state is advanced
    op   = @(x,mode) opSign_intrnl(m,n,seed,x,mode);

 case 2
    A  = (2.0 * (randn(m,n) < 0) - 1) / sqrt(m);
    op = opMatrix(A,{'Sign'});

 case 3
    seed = randn('state');
    for i=1:m, randn(n,1); end; % Ensure random state is advanced
    op = @(x,mode) opSignScaled_intrnl(m,n,seed,x,mode);

  otherwise
     error('Invalid mode')
end



function y = opSign_intrnl(m,n,seed,x,mode)
checkDimensions(m,n,x,mode);

% Store current random number generator state
seed0 = randn('state');
randn('state',seed);

if mode == 0
   y = {m,n,[0,1,0,1],{'Sign'}};
elseif mode == 1
   y = zeros(m,1);
   for i=1:n
      v = 2.0 * (randn(m,1) < 0) - 1;
      y = y + v * x(i);
   end
else
   y = zeros(n,1);
   for i=1:n
      v    = 2.0 * (randn(1,m) < 0) - 1;
      y(i) = v * x;
   end
end

% Restore original random number generator state
randn('state',seed0);



function y = opSignScaled_intrnl(m,n,seed,x,mode)
checkDimensions(m,n,x,mode);

% Store current random number generator state
seed0 = randn('state');
randn('state',seed);

if mode == 0
   y = {m,n,[0,1,0,1],{'Sign','Normalized'}};
elseif mode == 1
   y = zeros(m,1);
   for i=1:n
      v = 2.0 * (randn(m,1) < 0) - 1;
      y = y + v * x(i);
   end
   y = y / sqrt(m);
else
   y = zeros(n,1);
   for i=1:n
      v    = 2.0 * (randn(1,m) < 0) - 1;
      y(i) = v * x;
   end
   y = y / sqrt(m);
end

% Restore original random number generator state
randn('state',seed0);
