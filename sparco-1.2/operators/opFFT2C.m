function op = opFFT2C(m,n)
% OPFFT2C  Centralized two-dimensional fast Fourier transform (FFT).
%
%    OPFFT2C(M,N) creates a two-dimensional normalized Fourier
%    Transform that operates on matrices of size M by N. The
%    components are shifted to have to zero-frequency component in
%    the center of the spectrum. 

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: opFFT2C.m 1040 2008-06-26 20:29:02Z ewout78 $

op = @(x,mode) opFFT2C_intrnl(m,n,x,mode);


function y = opFFT2C_intrnl(m,n,x,mode)
checkDimensions(m*n,m*n,x,mode);
if (mode == 0)
   y = {m*n,m*n,[1,1,1,1],{'FFT2C'}};
elseif (mode == 1)
   y = fftshift(fft2(ifftshift(reshape(x,m,n)))) / sqrt(m*n);
   y = reshape(y,m*n,1);
else
   y = fftshift(ifft2(ifftshift(reshape(x,m,n)))) * sqrt(m*n);
   y = reshape(y,m*n,1);
end
