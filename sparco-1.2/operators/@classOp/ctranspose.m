function result = ctranspose(A)

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: ctranspose.m 1040 2008-06-26 20:29:02Z ewout78 $

A.adjoint = xor(A.adjoint,1);
result = A;

