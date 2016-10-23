function LockAxes(a)
% LockAxes -- Version-independent axis command
%  Usage
%    LockAxes(a)
%  Inputs
%    a     axis parameter, just as required by axis()
%
%  Side Effects
%    The axes are set to a and held, using a method
%    which works under both v3.5 and v4.0 of MATLAB.
%
% See Also
%    UnlockAxes, MATLABVERSION
%
global MATLABVERSION
if MATLABVERSION == 3.5,
   axis(a);
   plot(a(1)-.5, a(3)-.5);
   hold on;
else
   plot(a(1)-.5, a(3)-.5);
   axis(a);
   hold on;
end
       
    
%
% Copyright (c) 2006. David Donoho
%  

%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
