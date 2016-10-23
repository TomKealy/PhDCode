function [short_preamble,Xs]=short_train_seq(Nfft)
if nargin<1, Nfft=64; end
sqrt136 = sqrt(13/6)*(1+i);
Xs([3 11 23 39 43 47 51]) = sqrt136; % Frequency domain for k=-26:26
Xs([7 15 19 31 35]) = -sqrt136; 
xs = ifft([Xs(27:53) zeros(1,11) Xs(1:26)]); % Time domain
short_preamble = [xs xs xs(1:Nfft/2)];