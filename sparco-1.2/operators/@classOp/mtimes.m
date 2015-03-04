function y = mtimes(A,x)

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: mtimes.m 1040 2008-06-26 20:29:02Z ewout78 $

if isnumeric(A)
   y = (x' * A')';
else
   if A.adjoint == 0 % A*x
      if ismethod(A.op,'mtimes')
         y = A.op * x;
      else
         y = A.op(x,1);
      end

      if ~isempty(A.counter)
        evalin('base',[A.counter '=' A.counter '+ [1 0];']);
      end
   else
      if ismethod(A.op,'mtimes')
         if isnumeric(A)
           y = (x' * A.op)';
         else
           y = A.op' * x;
         end
      else
         y = A.op(x,2);
      end
   
      if ~isempty(A.counter)
         evalin('base',[A.counter '=' A.counter '+ [0 1];']);
      end
   end

   if A.sign < 0
      y = -y;
   end
end

