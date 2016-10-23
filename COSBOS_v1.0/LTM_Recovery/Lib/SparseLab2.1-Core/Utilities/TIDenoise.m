function y = TIDenoise(x,type,qmf)
% TIDenoise: Translation invariant denoising of a 1-D signal
% Usage
%   y = TIDenoise(x,type,qmf)
% Inputs
%   x	   1-d Noisy signal, dyadic
%   type   'S' for soft thresholding, 'H' for hard thresholding.
%		   Optional, default = hard thresholding
%   qmf    Quadrature mirror filter for wavelet transform.
%		   Optional, default = Symmlet 8
% Outputs
%   y      Denoised signal.
% References
%   Coifman, R.R. & Donoho, D.L. (1995). Translation-invariant
%     de-noising. In "Wavelets and Statistics", Antoniadis, A. &
%     Oppenheim, G. (Eds.), Lect. Notes Statist., 103, pp. 125-150,
%     New York: Springer-Verlag.

if nargin < 3,
    qmf = MakeONFilter('Symmlet',8);
end

if nargin < 2,
    type = 'H';
end

n = length(x);
lev = floor(log2(log(n)))+1;
thr = sqrt(2*log(n*log2(n)));

% Normalize noise level
[signalnorm, coef] = NormNoise(x,qmf);

% Compute translation invariant wavelet expansion
[tiwt] = FWT_TI(signalnorm, lev,qmf);
[nrow,ncol]  = size(tiwt);

if strcmp(type,'H'),
    tiwt(:,2:ncol) = HardThresh(tiwt(:,2:ncol),thr);
else
    tiwt(:,2:ncol) = SoftThresh(tiwt(:,2:ncol),thr);
end

y = (1/coef)*IWT_TI(tiwt,qmf);

    
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
