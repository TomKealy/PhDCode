%test_DSTBC_G2_PSK.m
clear, clf
MaxIter=1000;
sq2= sqrt(2); 
M=4; phs= 2*pi*[0:M-1]/M; 
MPSKs= exp(j*phs); % 4-PSK symbols
h= [0.9*exp(j*0.5) 1.1*exp(-j*0.4)]; % Assumed channel response
dh= 0.01*[exp(-j*0.1) -exp(j*0.2)]; % Channel variation
noise_amp= 0.2; % Amplitude of additive Gaussian noise 
x=[]; % Initialize the sequence of transmitted signal
S= eye(2); % S0=I;
dSh= 1;
r= h*S; % Eq.(9.4.67)
subplot(221)
for n=1:2:MaxIter
   xn= MPSKs(randint(1,2,M)+1); % +1,+j,-1, or -j
   % Constellation diagram of transmitted signal
   plot(real(xn),imag(xn),'ro'), hold on
   x= [x xn]; % Sequence of transmitted signal
   s= S*xn.'/sq2; % Encoded signal Eq.(9.4.66)
   S= [s [-s(2)'; s(1)']]; % Eq.(9.4.65)
   % Magnitude of transmitted signal 
   mag_of_transmitted_signal((n+1)/2)=sqrt(real(det(S))); 
   R= [r(1)' r(2); r(2)' -r(1)]; % Previously received signal
   detR= real(det(R));
   % Frequency-selective channel and additive Gaussian noise
   r= (h+dh*n)*S +noise_amp*randn(1,2); % Received signal Eq.(9.4.67)
   xhn= sq2/(-detR)*R*[r(1); r(2)']; % Decoding by Eq.(9.4.72)
   % Constellation diagram of decoded signal
   plot(xhn,'*'); hold on
   xhn= PSK_slicer(xhn,M); 
   xh([n n+1])= xhn; % Detection
end
ser= sum(abs(x-xh)>0.01)/MaxIter % Symbol error rate
plot(xh,'m+')
axis([-2 2 -2 2])
nn=0:length(mag_of_transmitted_signal)-1;
subplot(223)
plot(nn,abs(mag_of_transmitted_signal))