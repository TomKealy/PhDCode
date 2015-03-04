function f = pwsmoothfield(n,var,alpha)
% PWSMOOTHFIELD(N,VAR,ALPHA)
%   Generate an image of piecewise smooth filtered white noise
%
%   N = sidelength of the field in pixels
%   VAR = variance of original white nose
%   ALPHA = fraction of FFT coefficents to keep
%
%   Returns an N-by-N array
%
%   This file is used with the kind permission of Stephen J. Wright
%   (swright@cs.wisc.edu), and was originally included in the GPSR
%   v1.0 distribution: http://www.lx.it.pt/~mtf/GPSR .

% $Id: pwsmoothfield.m 1040 2008-06-26 20:29:02Z ewout78 $

f = sqrt(var)*randn(n);
F = fft2(f);
a = ceil(n*alpha/2);
b = fix(n*(1-alpha));
F(a+1:a+b,:) = 0;
F(:,a+1:a+b) = 0;
f = real(ifft2(F));

for i = 1:n
    for j = 1:n
        if (j/n >= 15*(i/n - 0.5)^3 + 0.4)
            f(i,j) = f(i,j) + sqrt(var);
        end
    end
end

