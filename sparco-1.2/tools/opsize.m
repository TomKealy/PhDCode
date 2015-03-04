function varargout = opsize(op,dim)
%OPSIZE  Get the dimensions of the operator
%
%   D = OPSIZE(OP) returns the dimensions of the operators as a row
%   vector.
%
%   [M,N] = OPSIZE(OP) returns the number of rows and columns of
%   the operator in M and N respectively.
%
%   OPSIZE(OP,DIM) returns the size of the operator in the
%   dimension DIM. For number of rows set DIM = 1, for number of
%   columsn set DIM = 2.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: opsize.m 1040 2008-06-26 20:29:02Z ewout78 $

% Check arguments
if nargin > 2
    error('Too many input arguments.')
end

info = op([],0);

if nargin == 2
    if ~isscalar(dim) || (dim < 1) || (round(dim) ~= dim)
        error('Dimension argument must be a positive integer scalar.');
    end
    
    if nargout > 1
        % One dimension queries, multiple outputs
        error('Unknown command option.');
    end
    
    if dim > 2
        varargout(1) = {1};
    else
        varargout(1) = {info{dim}};
    end
else
   if nargout <= 1
       varargout(1) = {[info{1},info{2}]};
   else
       varargout(1) = {info{1}};
       varargout(2) = {info{2}};
       for i=3:nargout, varargout(i) = {1}; end;
   end
end
