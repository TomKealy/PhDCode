function op = opBinary(m,n,mode)
%OPBINARY Binary (0/1) ensemble
%
%   OP = OPBINARY(M,N,MODE) creates an M by N binary ensemble
%   operator. When choosing MODE = 0 an explicit binary matrix is
%   formed and used when applying the operator. Choosing MODE = 1
%   causes the operator to generate the above matrix on the fly,
%   each time the operator is used. This mode can be used when
%   dealing with very large matrices.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: opBinary.m 1052 2008-07-02 20:52:08Z ewout78 $

if nargin < 3, mode = 0; end;


switch (mode)
 case 0
    A  = 1.0 * (randn(m,n) < 0);
    op = opMatrix(A,{'Binary'});

 case 1
    seed = randn('state');
    for i=1:m, randn(n,1); end; % Ensure random state is advanced
    op   = @(x,mode) opBinary_intrnl(m,n,seed,x,mode);

  otherwise
     error('Invalid mode')
end



function y = opBinary_intrnl(m,n,seed,x,mode)

% Store current random number generator state
seed0 = randn('state');
randn('state',seed);

checkDimensions(m,n,x,mode);

if mode == 0
   y = {m,n,[0,1,0,1],{'Binary'}};
elseif mode == 1
   y = zeros(m,1);
   for i=1:n
      v = 1.0 * (randn(m,1) < 0);
      y = y + v * x(i);
   end
else
   y = zeros(n,1);
   for i=1:n
      v    = 1.0 * (randn(1,m) < 0);
      y(i) = v * x;
   end
end

% Restore original random number generator state
randn('state',seed0);
