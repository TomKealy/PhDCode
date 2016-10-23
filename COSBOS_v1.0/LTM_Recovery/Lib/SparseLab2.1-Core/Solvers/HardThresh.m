function x = HardThresh(y,t)
% HardThresh -- Apply Hard Threshold 
%  Usage 
%    x = HardThresh(y,t)
%  Inputs 
%    y     Noisy Data 
%    t     Threshold
%  Outputs 
%    x     y 1_{|y|>t}
%
	x   = y .* (abs(y) > t);
  
%
% Copyright (c) 2006. Iddo Drori
%  

%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
