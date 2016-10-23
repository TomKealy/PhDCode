function y = SoftThresh(x,t)
s = abs(x) - t;
s = (s + abs(s))/2;
y = sign(x).*s;%
    
%
% Copyright (c) 2006. Iddo Drori
%  

% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
