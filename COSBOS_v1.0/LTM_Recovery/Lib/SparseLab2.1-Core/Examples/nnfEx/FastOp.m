function y = FastOp(mode, p, x)

global I;

if mode == 1
    
    z = fft(x);
    re = real(z);
    im = imag(z);
    y = [re(I); im(I)];
    
elseif mode == 2
    
    re = zeros(p,1);
    im = zeros(p,1);
    n = length(x);
    re(I) = x(1:n/2);
    im(I) = x(n/2+1:n);
    z = re + sqrt(-1)*im;
    y = real(ifft(z));
    
end

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
