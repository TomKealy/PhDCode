function op = opMask(A)
% OPMASK  Selection mask
%
%    OPMASK(A) creates an operator that computes the dot-product of
%    a given vector with the (binary) mask provided by A. If A is a
%    matrix it will be vectorized prior to use.
%
%    See also opColRestrict, opRestriction.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: opMask.m 1040 2008-06-26 20:29:02Z ewout78 $

op = @(x,mode) opMask_intrnl(A,x,mode);


function y = opMask_intrnl(A,x,mode)
m = size(A,1);
n = size(A,2);

checkDimensions(m*n,m*n,x,mode);
if mode == 0
   c =~isreal(A);
   y = {m*n,m*n,[c,1,c,1],{'Mask'}};
else
   y = A(:) .* x;
end
