function op = opHadamard(n,normalized)
% OPHADAMARD  Hadamard matrix
%
%    OPHADAMARD(N,NORMALIZED) creates a Hadamard operator for
%    vectors of length N. When N is not a power of two, vectors are
%    internally padded by zeroes to the next power of
%    two. Multiplication is done using a fast routine. When the
%    normalized flag is set, the columns of the Hadamard matrix are
%    scaled to unit two-norm. By default the NORMALIZED flag is set
%    to 0.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: opHadamard.m 1040 2008-06-26 20:29:02Z ewout78 $

if nargin < 2, normalized = 0; end;

op = @(x,mode) opHadamard_intrnl(n,normalized,x,mode);


function y = opHadamard_intrnl(n,normalized,x,mode)
checkDimensions(n,n,x,mode);
if mode == 0
   y = {n,n,[0,1,0,1],{'Hadamard',normalized}};
else
   s = length(x);
   y = opHadamardMult_intrnl([x;zeros(power(2,ceil(log2(s)))-s, 1)]);
   y = y(1:s);
   if normalized, y = y / sqrt(s); end;
end

function y = opHadamardMult_intrnl(x)
y = x;
n = length(y);
k = round(log2(n));
b = 1;     % Blocks on current level
s = n / 2; % Stride
for i=1:k  % Level
   for j=0:b-1  % Blocks
      for k=1:s   % Elements within block
         i1 = j*n + k;
         i2 = i1 + s;
         t1 = y(i1);
         t2 = y(i2);
         y(i1) = t1 + t2;
         y(i2) = t1 - t2;
      end
  end
  b = b * 2; s = s / 2; n = n / 2;
end
