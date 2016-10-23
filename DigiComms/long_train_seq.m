function [long_preamble,Xl]=long_train_seq(Nfft)
if nargin<1, Nfft=64; end
Xl([1:26 28:53]) = 1; % Frequency domain for k=-26:26
Xl([3 4 7 9 16 17 20 22 29 30 33 35 37:41 44 45 47 49])=-1;
xl = ifft([Xl(27:53) zeros(1,11) Xl(1:26)]); % Time domain
long_preamble= [xl(Nfft/2+1:Nfft) xl xl];