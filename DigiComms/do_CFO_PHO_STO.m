%do_CFO_PHO_STO.m
% To see the CFO/PHO/STO estimation and compensation
% Copyleft: Won Y. Yang, wyyang53@hanmail.net, CAU for academic use only
clear, clf
CFO=0.45; PHO=0.3;  Amp_noise=0.02; % CFO, PHO, and Noise amplitude
a = 10;  % Set a = 1~10 for full/no consideration of channel effect 
KC_CFO = 1; % Set to 1 or 0 for CFO compensation or not
KC_PHO = 1; % Set to 1 or 0 for Phase Offset compensation
Nfft=64; Ng=16; Ng1=Ng+1; Nfft1=Nfft+1; % N_Prefix=16; 
Nsym=Nfft+Ng; Nsym1=Nsym+1; Nw=Nsym; Nw1=Nw+1; Nw2=Nw*2; Nw4=Nw*4;  
Nbpsc=2; Ncp=4; Nnull=Nsym; N_OFDM=2;
ind_pilots=[8 22 44 58]; 
ind_excluding_pilots=[1:7 9:21 23:43 45:57 59:64];
gsymbols = 'v<>^'; %['s+xs'; 'd<^d'; '*v>*'; 's+xs']%N_OFDM symbols
ch_coef = exp(-a*[0:15]).'; H = fft(ch_coef,Nfft).';  % Channel response 
w_energy=zeros(1,Nw1); power=zeros(1,Nw1); energy_ratios=zeros(1,Nw1+31);
corr_sig1 = zeros(1,Ng1);   corr_sig2 = zeros(1,Nfft1);
ch_length = length(ch_coef);  ch_length1 = ch_length - 1;
ch_buf = zeros(1,ch_length);  rx_buf = zeros(1,Nsym1);  rxs_buf = [];
subplot(331), hold on
n = 0; nose = 0; nots = 0; % Numbers of symbol errors and total symbols
w_corr1 = 0; w_corr10 = 0; w_corr2 = 0; 
Max_energy_ratio=0; Min_energy_ratio=1e10; % 0.065
F_start = -Nsym; OFDM_start_points = [-Nsym];   
[short_preamble,S] = short_train_seq(Nfft);
[long_preamble,L] = long_train_seq(Nfft);
% Remaining part of previous frame, Null symbol, Time-domain preambles
Nd=80; any_symbol = 0.2*(rand(1,Nd)-0.5+j*(rand(1,Nd)-0.5));
tx = [any_symbol zeros(1,Nnull) short_preamble long_preamble];
Start_found = 0; Null_found = 0; Null_start_point = 0;
pm = 1;  pi2 = pi*2; % Pilot polarity 
coarse_CFO_est = 0;
for i=1:6+N_OFDM
   if i<7,  tx_cp = tx((i-1)*Nsym+[1:Nsym]);   
    else  msg = randint(1,Nfft*Nbpsc);
          if i>7,  tx_mod_prev = tx_mod;  
           elseif i==7,  tx_mod_prev = tx_cp; % long_tr not to be compared 
          end    % tx_mod transmitted previously              
          tx_mod=modulate_PSK_or_QAM(msg,Nbpsc,'QAM','bin'); % QAM mod
          tx_mod(ind_pilots) = [pm -pm pm pm];    % Pilot symbol
          tx_ifft = ifft(tx_mod,Nfft);              % IFFT
          tx_cp = [tx_ifft(end-Ng+1:end) tx_ifft];  % Add CP
          for m=ind_excluding_pilots
             txm = tx_mod(m); tmp = [real(txm) imag(txm)]>0;
             gsymbol = gsymbols(bin2deci(tmp,Nbpsc)+1); plot(txm,gsymbol)  
          end
   end
   for in=1:Nsym
      n = n+1;
      ch_buf = [tx_cp(in) ch_buf(1:end-1)];  % Channel buffer
      rxn=ch_buf*ch_coef*exp(j*(2*pi*CFO/Nfft*(n-1)+PHO)); % Add CFO & PHO
      rxn = rxn + Amp_noise*(randn+j*randn); rxs_buf = [rxs_buf rxn];
      rx_buf = [rx_buf(2:end) rxn]; % Received signal buffer
      power = [power(2:end) rx_buf(end)'*rx_buf(end)]; % Power window
      w_energy = [w_energy(2:end) w_energy(end)+power(end)]; % Energy window
      if n>Nw, w_energy(end) = w_energy(end)-power(1); end % Ng+Nfft
      if n>Nw1+30
        energy_ratio = w_energy(end)/w_energy(1);
        energy_ratios = [energy_ratios energy_ratio];
        if energy_ratio<Min_energy_ratio % 0.065?
          Min_energy_ratio = energy_ratio; 
         elseif Null_found<1, Null_start_point = n-1-Nw; Null_found = 1;
        end
        if energy_ratio>Max_energy_ratio
          Max_energy_ratio = energy_ratio; Start_found = 0; 
         elseif Start_found<10&Null_start_point>0&n-Null_start_point-Nw>50
          F_start = n-Start_found-Nw; Start_found = Start_found+1;
        end
      end
      corr_sig1(1:end-1)=corr_sig1(2:end); 
      corr_sig2(1:end-1)=corr_sig2(2:end);
      if F_start>0&n>F_start
        % Correlation between 2 points across Ng samples for short preamble
        if n<F_start+Nw2        
          corr_sig1(end) = rx_buf(end)*rx_buf(end-Ng)';  
          w_corr1 = w_corr1 + corr_sig1(end); 
        end  
        if n>F_start+Ng,  w_corr1 = w_corr1 - corr_sig1(1);  end
        if n==F_start+8*Ng, w_corr10 = w_corr1;  end
        if n==F_start+9*Ng
          coarse_CFO_est = (angle(w_corr10)+angle(w_corr1))/pi; 
        end
      end
      if Start_found>9&n>F_start+Nw2
        % Correlation between 2 points across Nfft samples for long preamble
        if KC_CFO>0
          rx_buf(end)=rx_buf(end)*exp(-j*2*pi*coarse_CFO_est/Nfft*(n-1)); 
        end
        if n<=F_start+Nw4
          corr_sig2(end) = rx_buf(end)*rx_buf(end-Nfft)';  
          w_corr2 = w_corr2 + corr_sig2(end); 
        end
        if n>F_start+Nw2+Nfft,  w_corr2 = w_corr2 - corr_sig2(1);  end
        if n==F_start+Nw4, fine_CFO_est = angle(w_corr2)/pi2;  end
      end
      if F_start>0&n>F_start+Nw4+Nw&(mod(n,Nsym)==2)  
        rx_wo_CP = rx_buf(end-Nfft1:end-2); % CP removed
        if KC_CFO>0  % CFO compensation
          rx_wo_CP=rx_wo_CP.*exp(-j*2*pi*fine_CFO_est/Nfft*[n-2-Nfft:n-3]); 
        end  
        rx_fft = fft(rx_wo_CP,Nfft);
        rx_fft_ch_comp = rx_fft./H;  % Channel compensation 
        ph_detected = angle(rx_fft_ch_comp(ind_pilots)*[pm;-pm;pm;pm]);
        rx_fft_ph_comp = rx_fft_ch_comp;
        if KC_PHO>0, rx_fft_ph_comp=rx_fft_ch_comp.*exp(-j*ph_detected); end 
        rx_sliced = QAM4_slicer(rx_fft_ph_comp);
        rx_dem=demodulate_PSK_or_QAM(rx_sliced,Nbpsc,'QAM','bin'); 
        tx_block1 = tx_mod_prev(ind_excluding_pilots);
        rx_sliced1 = rx_sliced(ind_excluding_pilots);
        nose = nose+sum(tx_block1~=rx_sliced1); nots = nots+length(tx_block1);
        for m= ind_excluding_pilots 
           txm = tx_mod_prev(m); tmp = [real(txm) imag(txm)]>0;
           gsymbol = gsymbols(bin2deci(tmp,Nbpsc)+1);
           subplot(332), plot(rx_fft_ch_comp(m),gsymbol), hold on
           subplot(333), plot(rx_sliced(m),gsymbol), hold on
        end
      end
   end
end
height=10;
subplot(312), plot(height*abs(rxs_buf),'b'), hold on, plot(energy_ratios,'r') 
True_start_points = 1+[0 80 160 320 480 560 640];
N_True_start_points=length(True_start_points);
stem(True_start_points,height*ones(1,N_True_start_points),'k*')
stem(F_start,1.2*height,'rx')
title('Estimated Frame Start Point and Energy Ratio across Nnull samples')
fprintf('\n Frame start %d has been detected to be %d.',Nd+Nsym1,F_start)
fprintf('\n The PHO %5.3f has been estimated to be %5.3f',PHO,ph_detected)
fprintf('\n The CFO %5.3f has been estimated to be', CFO)
if KC_CFO>0
  CFO_est = coarse_CFO_est + fine_CFO_est;   
  form=['\n %5.3f=%5.3f(coarse CFO estimate)+%5.3f(fine CFO estimate)\n'];
  fprintf(form,CFO_est,coarse_CFO_est,fine_CFO_est)
 else
  form=['\n %5.3f(coarse CFO estimate) and %5.3f(fine CFO estimate).\n'];
  fprintf(form,coarse_CFO_est,fine_CFO_est)
end
fprintf('The SER turned out to be %6.4f=%d/%d.\n',nose/nots,nose,nots)