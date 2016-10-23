%do_STO.m
% To see the effect of STO (symbol time offset) and its estimation/compensation
% Copyleft: Won Y. Yang, wyyang53@hanmail.net, CAU for academic use only
clear, clf
a = 10;  % Set a = 1~10 for full/no consideration of channel effect 
Amp_noise = 0.005; PHO = 0.3; % Noise amplitude and phase offset
Cor_thd = 0.924; Dif_thd = 0.065; % Threshold to detect the peak of correlation
KC_PH = 1;  % Set to 1 or 0 for Phase compensation or not
KC_STO = 1; % Set to 1/2/0 for ML/Classen STO estimation or exact symbol timing
if KC_STO==1, str1='Correlation';  str2='ML'; Thd=Cor_thd;
 elseif KC_STO==2, str1='Difference'; str2='Classen'; Thd=Dif_thd;
 else str1=''; str2='No STO estimation'; 
end  
Nfft=64; Ng=16; Ng1=Ng+1; Nfft1=Nfft+1; % N_Prefix=16; 
Nsym=Nfft+Ng; Nsym1=Nsym+1; Nnull=Nsym; Nnull1=Nnull+1; 
Nbpsc=2; Ncp=4; N_OFDM=5;   gsymbols = 'v<>^'; % N_OFDM symbols
True_starts=Ng1+[0:N_OFDM-1]*Nsym;
ch_coef = exp(-a*[0:15]).';  H = fft(ch_coef,Nfft).';  % Channel response 
w_energy = zeros(1,Nfft1);  power = zeros(1,Ng1);
corr_sig = zeros(1,Ng1);  corrs = zeros(1,Nsym);
ch_length = length(ch_coef); ch_length1 = ch_length - 1;
ch_buf = zeros(1,ch_length); rx_buf = zeros(1,Nsym1); rxs_buf = [];
pi2 = pi*2; pm = 1; % Pilot polarity
n = 0; i_STO = 0; nose = 0; win_corr = 0; % nose = Number of symbol errors
OFDM_starts = [-Nsym]; % An imagined previous OFDM symbol start point 
for i=1:N_OFDM
   msg = randint(1,Nfft*Nbpsc); % msgs = [msgs msg];
   if i>1,  tx_mod_prev = tx_mod;  end      % tx_mod transmitted previously
   tx_mod = modulate_PSK_or_QAM(msg,Nbpsc,'QAM','bin');
   tx_mod([8 22 44 58]) = [pm -pm pm pm];   % Pilot symbols
   tx_ifft = ifft(tx_mod,Nfft);             % IFFT
   tx_cp = [tx_ifft(end-Ng+1:end) tx_ifft]; % Add CP
   subplot(451), hold on, title('Transmitted symbols')
   for m=[1:7 9:21 23:43 45:57 59:64]
      txm = tx_mod(m); tmp = [real(txm) imag(txm)]>0;
      gsymbol = gsymbols(bin2deci(tmp,Nbpsc)+1);  plot(txm,gsymbol)
   end   
   for in=1:Nsym
      n = n+1;
      ch_buf = [tx_cp(in) ch_buf(1:end-1)];  % Channel buffer
      rxn = ch_buf*ch_coef*exp(j*PHO) + Amp_noise*(randn+j*randn); 
      rx_buf = [rx_buf(2:end) rxn]; rxs_buf = [rxs_buf rxn]; % Received signal
      % Window of instantaneous signal power
      power = [power(2:end) rx_buf(end)'*rx_buf(end)];  
      % Window of signal energy for duration Ng
      w_energy = [w_energy(2:end) w_energy(end)+power(end)]; 
      if n>Ng, w_energy(end) = w_energy(end)-power(end-Ng); end % Ng+Nfft
      corr_sig(1:end-1) = corr_sig(2:end);  
      if n>Nfft&KC_STO>0
        % Correlation/Difference between two points spaced Nfft samples apart
        if KC_STO==1
          corr_sig(end) = rx_buf(end)'*rx_buf(end-Nfft); % Correlation
         else 
          tmp = rx_buf(end)-rx_buf(end-Nfft); % Difference
          corr_sig(end) = tmp*tmp'; % Squared absolute difference
        end    
        win_corr = win_corr + corr_sig(end); 
      end
      % Windowed correlation for Ng points
      if n>Nsym, win_corr = win_corr-corr_sig(1); end
      if n>Nsym %+Ng  % CP-based Estimation of Symbol Start Time
        % Normalized, windowed correlation across Nfft samples for Ng points  
        corr_CP = abs(win_corr)/sqrt(w_energy(end)*w_energy(end-Nfft));  
        corrs = [corrs corr_CP];
        ellapsed = (n-Nsym+Ng > OFDM_starts(end)+Nfft);
        too_ellapsed = (n-Nsym+Ng > OFDM_starts(end)+Nsym+2);
        if KC_STO==1 % ML (Maximum Likelihood)
          condition1=(corr_CP>Cor_thd&corrs(n)<corrs(n-1)&ellapsed);
          detected = condition1|too_ellapsed; % Supposed maximum
         elseif KC_STO==2 % Classen
          condition2=(win_corr<Dif_thd&corrs(n)>corrs(n-1)&ellapsed);
          detected = condition2|too_ellapsed; % Supposed minimum
         else  
          detected = (mod(n,80)==2); % Exact symbol timing assumed
        end
        if detected  % If the detection has been done
          i_STO = i_STO + 1; ST_new=n-Nsym+Ng-(KC_STO~=2);
          OFDM_starts = [OFDM_starts ST_new];
          rx_wo_CP = rx_buf([end-Nfft1:end-2]+(KC_STO>0)); % CP removed
          rx_fft = fft(rx_wo_CP,Nfft);
          rx_fft_ch_comp = rx_fft./H;  % Channel compensation
          % Phase detection based on the pilot symboll
          ph_detected = angle(rx_fft_ch_comp([8 22 44 58])*[pm;-pm;pm;pm]);
          rx_fft_ph_comp = rx_fft_ch_comp;
          if KC_PH>0  % Phase compensation
            rx_fft_ph_comp = rx_fft_ch_comp.*exp(-j*ph_detected);
          end
          rx_sliced = QAM4_slicer(rx_fft_ph_comp);
          rx_dem = demodulate_PSK_or_QAM(rx_sliced,Nbpsc,'QAM','bin');
          % Extract the signals corresponding to the data subcarriers
          %  from the transmitted/received OFDM symbols for comparison
          tx_block1 = tx_mod_prev([1:7 9:21 23:43 45:57 59:64]);
          rx_sliced1 = rx_sliced([1:7 9:21 23:43 45:57 59:64]);
          nose1 = sum(tx_block1~=rx_sliced1); nose = nose + nose1;
          form1='\n # of symbol errors =%3d with STO estimate %3d for STO=%3d';
          fprintf(form1,nose1,ST_new,True_starts(i_STO))
          form2='\n The phase offset %5.3f has been estimated to be %5.3f.';
          fprintf(form2,PHO,ph_detected)
          if i_STO<5
            for m=[1:7 9:21 23:43 45:57 59:64] 
               txm = tx_mod_prev(m);
               tmp = [real(txm) imag(txm)]>0;
               gsymbol = gsymbols(bin2deci(tmp,Nbpsc)+1);    
               subplot(451+i_STO), plot(rx_fft_ch_comp(m),gsymbol), hold on
               if i_STO==3, title('Channel/STO-compensated signals'); end
               subplot(4,5,11+i_STO), plot(rx_sliced(m),gsymbol), hold on
               if i_STO==3, title('Constellation after slicing'); end
            end
          end
        end
     end
   end
end
no_of_symbol_errors = nose;
subplot(412), plot(corrs), hold on
plot(Thd*ones(size(corrs)),'r:') % plot the threshold line
fprintf('\nDetected start points: %d %d %d %d %d %d',OFDM_starts(2:end))
Detected_starts = OFDM_starts(2:end);
length_txs_cp = n; % Total number of samples
True_starts=[Ng1:Nsym:n];
fprintf('\nTrue starts: %d %d %d %d %d %d\n',True_starts)
stem(True_starts,1.2*ones(size(True_starts)),'k*')
N_Sym=length(Detected_starts);
stem(Detected_starts,1.3*ones(1,N_Sym),'rx')
title('Detected Start Points of Symbols and Correlation across Nfft samples')