function y = UpDyadLo(x,qmf)
% UpDyadLo -- Lo-Pass Upsampling operator; periodized
%  Usage
%    u = UpDyadLo(d,f)
%  Inputs
%    d    1-d signal at coarser scale
%    f    filter
%  Outputs
%    u    1-d signal at finer scale
%
%  See Also
%    DownDyadLo, DownDyadHi, UpDyadHi, IWT_PO, iconv
%
	y =  iconv(qmf, UpSample(x) );

    
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
