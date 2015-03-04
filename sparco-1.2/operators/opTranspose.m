function op = opTranspose(op)
% OPTRANSPOSE   Transpose operator
%
%    OPTRANSPOSE(OP) creates an operator that is the the adjoint
%    operator OP.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: opTranspose.m 1040 2008-06-26 20:29:02Z ewout78 $

op = @(x,mode) opTranspose_intrnl(op,x,mode);


function y = opTranspose_intrnl(op,x,mode)
info = op([],0);
checkDimensions(info{2},info{1},x,mode);
if mode == 0
  y = {info{2},info{1},info{3}([3,4,1,2]),{'Transpose',op}};
elseif mode == 1
  y = op(x,2);
else
  y = op(x,1);
end
