function UnlockAxes
% UnlockAxes -- Version-independent axis command
%  Usage
%    UnlockAxes
%
%  Side Effects
%    Cancels the *hold* side effect of LockAxes using a method
%    which works under both v3.5 and v4.0 of MATLAB.
%
%  See Also
%    LockAxes, MATLABVERSION
%
global MATLABVERSION
if MATLABVERSION == 3.5,
   hold off; axis;
else
   hold off;
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
