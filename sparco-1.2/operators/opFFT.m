function op = opFFT(n)
% OPFFT  One-dimensional fast Fourier transform (FFT).
%
%    OPFFT(N) create a one-dimensional normalized Fourier transform
%    operator for vectors of length N.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: opFFT.m 1040 2008-06-26 20:29:02Z ewout78 $

op = @(x,mode) opFFT_intrnl(n,x,mode);


function y = opFFT_intrnl(n,x,mode)
checkDimensions(n,n,x,mode);
if mode == 0
   y = {n,n,[1,1,1,1],{'FFT'}};
elseif mode == 1
   y = fft(x) / sqrt(length(x));
else
   y = ifft(x) * sqrt(length(x));
end
