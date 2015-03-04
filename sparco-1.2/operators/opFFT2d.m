function op = opFFT2d(m,n)
% OPFFT2D  Two-dimensional fast Fourier transform (FFT)
%
%    OPFFT2D(M,N) creates a two-dimensional normalized Fourier
%    transform operator for matrices of size M by N. Input and
%    output of the matrices is done in vectorized form.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: opFFT2d.m 1040 2008-06-26 20:29:02Z ewout78 $

op = @(x,mode) opFFT2d_intrnl(m,n,x,mode);


function y = opFFT2d_intrnl(m,n,x,mode)
checkDimensions(m*n,m*n,x,mode);
if mode == 0
   y = {m*n,m*n,[1,1,1,1],{'FFT2d'}};
elseif mode == 1
   y = reshape( fft2(reshape(x,m,n)) / sqrt(m*n), m*n, 1);
else
   y = reshape(ifft2(reshape(x,m,n)) * sqrt(m*n), m*n, 1);
end
