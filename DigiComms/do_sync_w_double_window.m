%do_sync_w_double_window.m
clear, clf
Cor_thd=0.988; % The threshold to determine the peak of correlation
Nfft=64; Ng=16; Nsym=Nfft+Ng; Nsym1=Nsym+1; 
Nnull=Nsym; Nw=Nnull; Nw1=Nw+1; Nw2=Nw*2; Ng1=Ng+1; Nfft1=Nfft+1; 
Nd=90; % Remaining period of the last symbol in the previous frame
N_OFDM=3; % One Null + N_OFDM symbols
Max_energy_ratio=0;  Min_energy_ratio=1e10;
r = [rand(1,Nd)-0.5+j*(rand(1,Nd)-0.5) zeros(1,Nnull)];
for i=1:N_OFDM
   symbol=rand(1,Nfft)-0.5 +j*(rand(1,Nfft)-0.5); 
   r = [r symbol(end-Ng+1:end) symbol];
end
L_frame = length(r); 
r = r + 0.1*(rand(1,L_frame)-0.5+j*(rand(1,L_frame)-0.5));
energy_w1=zeros(1,Nw1); power_w=zeros(1,Nw1); 
sig_w=zeros(1,Nfft1); energy_w2=zeros(1,Nfft1); corr_w=zeros(1,Ng1);
OFDM_start_points = [0]; corr=0; 
for n=1:L_frame
   sig_w = [sig_w(2:end) r(n)]; % Signal window
   power_n = r(n)'*r(n); % Current signal power
   power_w = [power_w(2:end)  power_n]; % Power window 
   energy_w1=[energy_w1(2:end) energy_w1(end)+power_n]; %Energy window
   if n>Nw, energy_w1(end)=energy_w1(end)-power_w(1); end %of size Nw
   energy_w2=[energy_w2(2:end) energy_w2(end)+power_n]; %Energy window
   if n>Ng, energy_w2(end)=energy_w2(end)-power_w(end-Ng); end
   corr_w(1:end-1) = corr_w(2:end); 
   if n>Nfft
     %Correlation between signals at 2 points spaced Nfft samples apart
     corr_w(end)=abs(sig_w(end)'*sig_w(1)); corr=corr+corr_w(end);
   end
   if n>Nsym, corr=corr-corr_w(1); end %Correlation window of Ng pts
   % Null Symbol detection based on energy ratio
   if n>=Nw2
     energy_ratio = energy_w1(end)/energy_w1(1);
     energy_ratios(n) = energy_ratio; 
     if energy_ratio<Min_energy_ratio  % Eq.(11.3.8a)
       Min_energy_ratio = energy_ratio; Null_start_point = n-Nw+1;
     end
     if energy_ratio>Max_energy_ratio  % Eq.(11.3.8b)
       Max_energy_ratio = energy_ratio; F_start = n-Nw+1;
     end
   end
   % CP-based Symbol Time estimation
   if n>Nsym
     % Normalized, windowed correlation across Nfft samples for Ng pts  
     correlation=corr/sqrt(energy_w2(end)*energy_w2(1)); % Eq.(11.3.9)
     correlations(n) = correlation; 
     if correlation>Cor_thd&n-Nsym>OFDM_start_points(end)+Nfft
       OFDM_start_points = [OFDM_start_points n-Nsym+1];
     end
   end
end
Estimated_start_points=[Null_start_point F_start OFDM_start_points(2:end)]
True_start_points=[Nd+1:Nsym:L_frame]
N_True_start_points=length(True_start_points);
subplot(311), stem(real(r)), set(gca,'XTick',True_start_points)
hold on, stem(True_start_points,0.9*ones(1,N_True_start_points),'k*')
N_Sym=length(Estimated_start_points);
stem(Estimated_start_points,1.1*ones(1,N_Sym),'rx')
title('Estimated Start Points of Symbols')
subplot(312),semilogy(energy_ratios),set(gca,'XTick',True_start_points)
title('Ratio of 2 Successive Windowed Energies for Nw samples')
subplot(313), plot(correlations)
hold on, title('Correlation across Nfft samples')