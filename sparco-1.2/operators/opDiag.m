function op = opDiag(n,d)
% OPDIAG   Scaling operator (i.e., diagonal matrix).
%
%    OPDIAG(N,D) creates an operator for N by N diagonal matrix
%    multiplication. The diagonal is given by D which can be either
%    a vector of length N or a scalar. OPDIAG(D) is a shortcut for
%    OPDIAG(LENGTH(D),D) and can be used instead.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: opDiag.m 1040 2008-06-26 20:29:02Z ewout78 $

if nargin < 2
  d = n;
  n = length(d);
end

d  = d(:); % Ensure d is a column vector
op = @(x,mode) opDiag_intrnl(n,d,x,mode);


function y = opDiag_intrnl(n,d,x,mode)
checkDimensions(n,n,x,mode);
if (mode == 0)
   y = {n,n,[0,1,0,1],{'Diag'}};
elseif (mode == 1)
   y = d.*x;
else
   y = conj(d).*x;
end
