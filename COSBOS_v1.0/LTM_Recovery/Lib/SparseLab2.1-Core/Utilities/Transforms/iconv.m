function y = iconv(f,x)
% iconv -- Convolution Tool for Two-Scale Transform
%  Usage
%    y = iconv(f,x)
%  Inputs
%    f   filter
%    x   1-d signal
%  Outputs
%    y   filtered result
%
%  Description
%    Filtering by periodic convolution of x with f
%
%  See Also
%    aconv, UpDyadHi, UpDyadLo, DownDyadHi, DownDyadLo
%
	n = length(x);
	p = length(f);
	if p <= n,
	   xpadded = [x((n+1-p):n) x];
	else
	   z = zeros(1,p);
	   for i=1:p,
		   imod = 1 + rem(p*n -p + i-1,n);
		   z(i) = x(imod);
	   end
	   xpadded = [z x];
	end
	ypadded = filter(f,1,xpadded);
	y = ypadded((p+1):(n+p));

%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
