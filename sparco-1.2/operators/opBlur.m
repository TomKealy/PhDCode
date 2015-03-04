function op = opBlur(m,n)
%OPBLUR   Two-dimensional blurring operator
%
%   OPBLUR(M,N) creates an blurring operator for M by N
%   images. This function is used for the GPSR-based test problems
%   and is based on the implementation by Figueiredo, Nowak and
%   Wright, 2007.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: opBlur.m 1040 2008-06-26 20:29:02Z ewout78 $

yc = round(m/2 + 1);
xc = round(n/2 + 1);

% Create blurring mask
h  = zeros(m,n);
for i=-4:4
   h(i+yc,(-4:4)+xc)= (1./(1+i*i+(-4:4).^2)); 
end
h = fftshift(h);   
h = h/sum(h(:));

op = @(x,mode) opBlur_intrnl(m,n,h,x,mode);


function y = opBlur_intrnl(m,n,h,x,mode)
checkDimensions(m*n,m*n,x,mode);
if mode == 0
   y = {m*n,m*n,[0,1,0,1],{'Blur'}};
elseif mode == 1
   y = reshape(ifft2(fft2(h).*fft2(reshape(x,m,n))),m*n,1);
   if isreal(x), y = real(y); end;
else
   y = reshape(ifft2(conj(fft2(h)).*fft2(reshape(x,m,n))),m*n,1);
   if isreal(x), y = real(y); end;
end
