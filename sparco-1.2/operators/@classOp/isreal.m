function r = isreal(A)

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: isreal.m 1040 2008-06-26 20:29:02Z ewout78 $

% Get size information
if ismethod(A.op,'isreal') || isnumeric(A.op)
   r = isreal(A.op);
else
   info = A.op([],0); c = info{3};
   if (c(1)==0) && (c(1)==0)
      r = 1;
   else
      r = 0;
   end
end
