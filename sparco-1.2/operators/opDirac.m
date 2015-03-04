function op = opDirac(n)
% OPDIRAC  Identity operator
%
%    OPDIRAC(N) creates the identity operator for vectors of length
%    N. 

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: opDirac.m 1040 2008-06-26 20:29:02Z ewout78 $

op = @(x,mode) opDirac_intrnl(n,x,mode);


function y = opDirac_intrnl(n,x,mode)
checkDimensions(n,n,x,mode);
if (mode == 0)
   y = {n,n,[0,1,0,1],{'Dirac'}};
else
   y = x;
end

