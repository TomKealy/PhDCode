function [s,varargout] = size(A,n)

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: size.m 1040 2008-06-26 20:29:02Z ewout78 $

% Check validity of parameter 'n'
if nargin > 1
  if ~isscalar(n) || ~((n == 1) || (n == 2))
    error('Dimension argument must be a positive integer scalar within indexing range.')
  end
end

% Get size information
if ismethod(A.op,'size') || isnumeric(A.op)
  opSize = size(A.op);
else
  opInfo = A.op([],0);
  opSize = [opInfo{1}, opInfo{2}];
end
if A.adjoint, opSize = opSize([2,1]); end

% Get dimension to return
if nargin > 1
  if n <= 2, opSize = opSize(n); else opSize = 1; end
end

if nargout <= 1
  s = opSize;
elseif nargout >= 2
  if nargin == 1
     s = opSize(1);
     varargout{1} = opSize(2);
     for i=2:nargout, varargout{i} = 1; end
  else
     error('Unknown command option.');
  end
end
