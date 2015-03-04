function r = opisreal(op)
%OPISREAL  True for real operator
%
%   OPISREAL(OP) returns 1 if the result of multiplication of real
%   vectors by OP and OP transpose is real, and 0 otherwise.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: opisreal.m 1040 2008-06-26 20:29:02Z ewout78 $

% Check arguments
if nargin > 1
    error('Too many input arguments.')
end

try
  info = op([],0); c = info{3};
  if (c(1)==0) && (c(1)==0)
     r = 1;
  else
     r = 0;
  end
catch
  error('Parameter must be an operator');
end
