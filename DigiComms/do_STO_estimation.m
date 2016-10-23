%do_STO_estimation.m 
% To estimate the STO (Symbol Time Offset)
%Copyleft: Won Y. Yang, wyyang53@hanmail.net, CAU for academic use only
clear, clf
Cor_thd = 0.93; % Threshold of correlation for peak detection
Nfft=64;  Ng=16; % FFT/Guard interval size
Nsym=Nfft+Ng;  Nw=Ng; % Symbol/Window size
Nd=65; % Remaining period of the last symbol in the previous frame
% Make a pseudo frame to transmit
t_frame = rand(1,Nd)-0.5+j*( rand(1,Nd)-0.5);  
N_Symbols = 3; % Number of OFDM symbols to be generated for simulation
for i=1:N_Symbols
   symbol = rand(1,Nfft)-0.5+j*(rand(1,Nfft)-0.5); % An arbitrary data
   symbol_cp = [symbol(end-Ng+1:end) symbol]; % OFDM Symbol with CP
   t_frame = [t_frame symbol_cp]; % Append a frame by a symbol with CP
end
L_frame = length(t_frame);  nn=0:L_frame-1;
noise = 0.2*(rand(1,L_frame)-0.5 +j*(rand(1,L_frame)-0.5));
r = t_frame + noise;
sig_w = zeros(2,Nw); % Initialize the two sliding window buffers
STOs = [0]; % Initialize the STO buffer
for n=1:L_frame
   sig_w(1,:) = [sig_w(1,2:end) r(n)]; % Update signal window 1
   m = n-Nfft;
   if m>0
     sig_w(2,:)=[sig_w(2,2:end) r(m)]; % Update signal window 2
     den = norm(sig_w(1,:))*norm(sig_w(2,:));
     corr(n) = abs(sig_w(1,:)*sig_w(2,:)')/den;
     if corr(n)>Cor_thd & m>STOs(end)+Nsym-15
       STOs=[STOs  m]; % List the estimated STO
     end
   end
end
Estimated_STOs = STOs(2:end)
True_STOs = Nd+Ng + [0:N_Symbols-1]*Nsym
subplot(311)
stem(nn,real(r)), ylim([-0.6  1.1])
hold on, stem(True_STOs,0.8*ones(size(True_STOs)),'k*')
stem(Estimated_STOs,0.6*ones(size(Estimated_STOs)),'rx')
title('Estimated Starting Times of OFDM Symbols')
subplot(312)
plot(nn+1,corr), ylim([0  1.2]), hold on, stem(Estimated_STOs+Nfft,corr(Estimated_STOs+Nfft),'r:^')
% The points at which the correlation is presumably maximized,
%  yielding the STO estimates.  
stem(Estimated_STOs,0.9*ones(size(Estimated_STOs)),'rx')
title('Correlation between two sliding windows across Nfft samples')
set(gca,'XTick',sort([Estimated_STOs Estimated_STOs+Nfft]))