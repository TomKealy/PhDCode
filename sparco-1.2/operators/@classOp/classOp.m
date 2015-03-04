function result = classOp(op,counter_name)
%CLASSOP   Turn operator into class
%   RESULT = CLASSOP(OP,COUNTER_NAME), creates a class instance
%   for operator OP and returns it in RESULT. Optional parameter
%   COUNTER_NAME can be used to specify which global variable to
%   use for counting the number of direct and adjoint operator
%   applications done.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: classOp.m 1040 2008-06-26 20:29:02Z ewout78 $

if nargin < 2
  counter_name = [];
end

result.adjoint =  0;
result.sign    = +1;
result.op      = op;
result.counter = counter_name;

if ~isempty(counter_name)
  evalin('base',[counter_name ' = [0,0];']);
end

result = class(result,'classOp');
