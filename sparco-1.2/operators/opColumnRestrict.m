function op = opColumnRestrict(m,n,idx,type)
% OPCOLUMNRESTRICT   Restriction operator on matrix columns
%
%    OPCOLUMNRESTRICT(M,N,IDX,TYPE), with TYPE = 'DISCARD', creates
%    an operator that extracts the columns indicated by IDX from a 
%    given M by N matrix. The adjoint operator takes an M by
%    length(IDX) matrix and outputs an M by N matrix with the
%    columns filled by the given input matrix. Note that all input
%    and output matrices are in vectorized from. When TYPE =
%    'ZERO', all columns that are not in IDX are zero-padded
%    instead of discarded.
%
%    See also opMask, opRestriction.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: opColumnRestrict.m 1040 2008-06-26 20:29:02Z ewout78 $

if (min(idx) < 1) || (max(idx) > n)
   error('Index parameter must be integer and match dimensions of the data');
end

if (nargin < 4) || isempty(type)
  type = 'discard';
end

if strcmp(lower(type),'discard')
   l  = length(idx);
   op = @(x,mode) opColumnRestrict_intrnl1(idx,l,m,n,x,mode);
else
   invidx = ones(n,1);
   invidx(idx) = 0;
   idx = find(invidx);
   op  = @(x,mode) opColumnRestrict_intrnl2(idx,m,n,x,mode);
end 


function y = opColumnRestrict_intrnl1(idx,l,m,n,x,mode)
checkDimensions(m*l,m*n,x,mode);
if mode == 0
   y = {m*l,m*n,[0,1,0,1],{'ColumnRestrict','discard'}};
elseif mode == 1
   y = reshape(x,m,n);
   y = reshape(y(:,idx),m*l,1);
else
   y = zeros(m,n);
   y(:,idx) = reshape(x,m,l);
   y = reshape(y,m*n,1);
end

function y = opColumnRestrict_intrnl2(invidx,m,n,x,mode)
checkDimensions(m*n,m*n,x,mode);
if mode == 0
   y = {m*n,m*n,[0,1,0,1],{'ColumnRestrict','zero'}};
else
   y = reshape(x,m,n);
   y(:,invidx) = 0;
   y = y(:);
end
