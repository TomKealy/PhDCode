function y = MirrorFilt(x)
% MirrorFilt -- Apply (-1)^t modulation
%  Usage
%    h = MirrorFilt(l)
%  Inputs
%    l   1-d signal
%  Outputs
%    h   1-d signal with DC frequency content shifted
%        to Nyquist frequency
%
%  Description
%    h(t) = (-1)^(t-1)  * x(t),  1 <= t <= length(x)
%
%  See Also
%    DyadDownHi
%

	y = -( (-1).^(1:length(x)) ).*x;

    
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
