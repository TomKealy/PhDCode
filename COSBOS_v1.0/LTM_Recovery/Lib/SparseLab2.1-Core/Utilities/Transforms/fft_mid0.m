function y=fft_mid0(x)
% fftmid0 -- 1d fft with argument [-pi,pi] (instead of [0,2pi])
% Usage:
%   Y = fft_mid0(X)
% Inputs:
%   X	Array(n) 
% Outputs:
%   Y   Array(n)
% Description:
%  Performs 1d fft with grid (-n/2):(n/2-1) on both time
%  and frequency side. 
%    y(k) = sum_{t=-n/2}^{n/2-1} exp(i 2pi/n kt) x(t) , (-n/2) <= k < n/2
%

y = fftshift(fft(fftshift(x)));

%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
