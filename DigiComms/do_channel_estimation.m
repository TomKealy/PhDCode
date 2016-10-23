%do_channel_estimation.m
clear, clf
% Read the time-domain channel response stored in a 64x2 matrix
load  ch_complex.dat; 
h=(ch_complex(:,1)+j*ch_complex(:,2)).'; % Time-domain channel response
Nfft=64; Ng=16; Nsym=Nfft+Ng; Nnull=Nsym; Nw=2*Ng; Tsym=4e-6;
A_sig=0.5; A_noise=0.01; % Amplitudes of signal and noise
Nd=120; % Remaining period of the last symbol in the previous frame
[s_preamble,Short] = short_train_seq(Nfft);
[l_preamble,Long] = long_train_seq(Nfft);
% Make a pseudo received sequence
t_frame=[A_sig*(rand(1,Nd)-0.5) zeros(1,Nnull) s_preamble l_preamble];
symbol=A_sig*(rand(1,Nd)-0.5); symbol=[symbol(end-Ng+1:end) symbol];
t_frame = [t_frame symbol];  L_frame = length(t_frame); 
noise = A_noise*(rand(1,L_frame)-0.5 +j*(rand(1,L_frame)-0.5)); 
r_frame = channel(t_frame,h) + noise;
True_STO = Nd+Nnull+Nsym*2+1 % Start point of the long preamble 
STO = True_STO
% STO estimation is critical to the performance of channel estimation
% To realize this, change the above line into STO=True_STO+1 or -1
H_est = channel_estimate(Long,r_frame(STO+[0:159]));
% The true frequency-domain channel response is obtained from the FFT 
%  of the time-domain channel (impulse) response.
H = fft(h,Nfft); % k=0(1):26(27)  27(28):37(38)  38(39):63(64) 
H_true = H([39:64 2:27]); % Arranged in -/+ frequency k=-26:-1  1:26 
discrepancy_H_est_and_H_true = norm(H_est-H_true)/norm(H_true)
% Let's see how the channel equalizer with the estimated channel 
%  response (H_est) works.
X = rand(1,52)-0.5+j*(rand(1,52)-0.5); % for k=[-26:-1  1:26]
X_arranged = [0 X(27:end) zeros(1,11) X(1:26)]; % in +/- frequency
x = ifft(X_arranged); x_CP = [x(49:64) x];  % IFFT and add CP
y = channel(x_CP,h); % Channel output to an arbitrary input with CP
Y = fft(y(17:80)); % Remove CP and FFT
Y=Y([39:64 2:27]); % Arranged in -/+ frequency for k=[-26:-1 1:26]
Yeq = equalizer_in_freq(Y,H_est);
discrepancy_X_and_Y = norm(X-Y)/norm(X) % With no channel compensation 
discrepancy_X_and_Yeq = norm(X-Yeq)/norm(X) % With channel equalizer