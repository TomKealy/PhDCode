function op = opSplitComplex(n)
% OPSPLITCOMPLEX  Complex to real and imaginary
%
%    OPSPLITCOMPLEX(N) creates a nonlinear operator that splits a
%    complex vector of length N into its real and imaginary parts
%    [real(x); imag(x)]. In transpose mode it combines the real and
%    imaginary parts into a complex vector.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: opSplitComplex.m 1040 2008-06-26 20:29:02Z ewout78 $

op = @(x,mode) opComplexSplit_intrnl(n,x,mode);

function y = opComplexSplit_intrnl(n,x,mode)
checkDimensions(2*n,n,x,mode);
if mode == 0
   y = {2*n,n,[0,0,1,1],{'SplitComplex'}};
elseif mode == 1
   y = [real(x); imag(x)];
else
   y = x(1:n) + sqrt(-1) * x(n+(1:n));
end
