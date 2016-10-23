function [y,coef] = NormNoise(x,qmf)
% NormNoise -- Estimates noise level, Normalize signal to noise level 1
%  Usage
%    [y,coef] = NormNoise(x,qmf)
%  Inputs
%    x     1-d signal
%    qmf   quadrature mirror filter
%  Outputs
%    y     1-d signal, scaled so wavelet coefficients
%          at finest level have median absolute deviation 1.
%    coef  estimation of 1/sigma
%
%  Description
%    This is required pre-processing to use any of the DeNoising
%    tools on naturally-occurring data.
%
%  See Also
%    WaveShrink, CPDeNoise,WPDeNoise
%
	u = DownDyadHi(x,qmf);
	s = median(abs(u));
	if s ~= 0
		y =  .6745 .* x ./s;
                coef = .6745 /s;
	else
		y = x;
                coef = 1;
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
