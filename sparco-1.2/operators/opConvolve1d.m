function op = opConvolve1d(kernel,offset,circular)
% OPCONVOLVE1D   One-dimensional convolution operator
%
%    OPCONVOLVE1D(KERNEL,OFFSET,CIRCULAR) creates a one-dimensional
%    convolution operator with the KERNEL function. The length of
%    the KERNEL vector must coincide with the signal length. OFFSET
%    gives the index in the KERNEL vector to be taken as its
%    center and is set to 1 by default (i.e. no shifting). The
%    CIRCULAR flag indicates whether the signal is periodic. If set
%    to any non-zero value convolution is wrapped around the vector
%    endpoints. By default CIRCULAR is set to 1.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: opConvolve1d.m 1040 2008-06-26 20:29:02Z ewout78 $

if nargin < 2, offset = 1;   end;
if nargin < 3, circular = 1; end;

% Ensure kernel is a column vector
kernel = kernel(:);

% Get kernel length
n = length(kernel);

% Shift and extend kernel
if circular
  kernel = [kernel(offset:end); kernel(1:offset-1)];
else
  kernel = [zeros(n-(offset-1),1); kernel; zeros(offset-1,1)];
end

% Precompute kernel in Fourier mode
kf = fft(kernel);

op = @(x,mode) opConvolve1d_intrnl(kernel,kf,circular,n,x,mode);



function y = opConvolve1d_intrnl(kernel,kf,circular,n,x,mode)
checkDimensions(n,n,x,mode);
if mode == 0
  y = {n,n,[~isreal(kernel),1,~isreal(kernel),1],{'Convolve1d'}};
elseif mode == 1
  if circular
    y = ifft(kf.*fft(x));
  else
    xf = fft([zeros(n,1); x]);
    y  = ifft(kf.*xf);
    y  = y(1:n);
  end
  
  if (isreal(x) && isreal(kernel)), y = real(y); end;
else
  if circular
    y = ifft(conj(kf).*fft(x));
  else
    x = [x; zeros(n,1)];
    y = ifft(conj(kf).*fft(x));
    y = y(n+1:end);
  end  

  if (isreal(x) && isreal(kernel)), y = real(y); end;
end
