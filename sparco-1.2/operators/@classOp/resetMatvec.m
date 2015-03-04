function resetMatvec(A)

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: resetMatvec.m 1040 2008-06-26 20:29:02Z ewout78 $

if ~isempty(A.counter)
   evalin('base',[A.counter '= [0 0];']);
end
