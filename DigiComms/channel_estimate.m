function H_est=channel_estimate(X,y)
% X = Known frequency-domain training symbol
% y = Time-domain output of the channel 
% H_est = Estimate of the channel response
if length(X)>52,  X = X([1:26 28:53]);  end   % for k=[-26:-1  1:26]
y1 = y(32+[1:64]);  y2 = y(96+[1:64]);
Y1 = fft(y1); Y1 = Y1([39:64 2:27]); % Arranged in the -/+ frequency 
Y2 = fft(y2); Y2 = Y2([39:64 2:27]); % Arranged in the -/+ frequency 
H_est = X.*(Y1+Y2)/2; % Eq.(11.4.3)