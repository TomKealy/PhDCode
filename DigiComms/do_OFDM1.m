%do_OFDM1.m
% No initialization (Neither channel nor training for channel estimation)
Ng=16; Nrow=16; % Length of prefix & Number of rows for block interleaving
RATE=6; [CR,DR,Nbpsc,Ncbps,Ndbps,pun_vec,Mod,pps]=OFDM_parameters(RATE);
b=Nbpsc; M=2^b; % Number of vits per subcarriers and Modulation order
Gc=[171 133]; K=size(Gc,1); % Generator matrix and # of encoder input bits
Gc_m=max(Gc.'); % Constraint length vector 
for i=1:length(Gc_m), Lc(i)=length(deci2bin1(oct2dec(Gc_m(i))));  end
trel=poly2trellis(Lc,Gc); Tbdepth=sum(Lc)*5; delay=Tbdepth*K;
N_Symbol=500; % Number of OFDM symbols per iteration of simulation
EbN0dBs=[7  10]; MaxIters=[100  1000]; N_EbN0dBs=length(EbN0dBs);
SigPower= 0.013; %0.0130386; % Normalized average signal power
EbN0s=10.^(EbN0dBs/10)*str2num(CR);
Amps_noise=sqrt(SigPower/2./(b*EbN0s)/2); % sqrt(N0/2) for I/Q components
EbN0dBs_t=0:0.1:15; EbN0s_t=10.^(EbN0dBs_t/10);  Target_no_of_error=100;
for i_EbN0=1:N_EbN0dBs
   nobe=0; notb=0; % Initialize Numbers of bit errors and total bits 
   for iter=1:MaxIters(i_EbN0)
      msgs=[];   decodeds=[];
      state=[]; m=[]; s=[]; in=[]; N_data=Ndbps/b;
      for i_sym=1:N_Symbol
         msg = deci2bin(randint(1,N_data,M),b);
         if i_sym==N_Symbol,  msg(end-20:end)=0;  end
         msgs= [msgs  msg];
         if RATE==0, coded_msg = msg; % No channel encoding
          else  [coded_msg,state] = convenc(msg,trel,state);
         end
         punctured_msg = puncture(coded_msg,pun_vec); 
         interleaved_msg = interleaving(punctured_msg,Nrow,Ncbps,Nbpsc); 
         modulated = modulate_PSK_or_QAM(interleaved_msg,b,Mod,'bin');
         pn = pps(mod(i_sym-1,127)+1); % Pilot polarity
         tx = ifft(add_vc(modulated,pn)); % Add virtual carriers and IFFT
         ch_in = [tx(end-Ng+1:end) tx]; % Add the prefix
         noise = randn(size(ch_in)) + j*randn(size(ch_in));
         ch_out = ch_out + Amps_noise(i_EbN0)*noise;
         fft_out = fft(ch_out(Ng+1:end)); % Remove the prefix and FFT
         sym_vc_removed = remove_vc(fft_out); % Remove the virtual carriers
         demodulated = demodulate_PSK_or_QAM(sym_vc_removed,b,Mod);
         deinterleaved = deinterleaving(demodulated,Nrow,Ncbps,Nbpsc); 
         deinterleaved(find(deinterleaved==0)) = -1; 
         deinterleaved = -deinterleaved; % Convert into negative bipolar
         depunctured = depuncture(deinterleaved,pun_vec); 
         if RATE==0,  decoded = depunctured; % No channel decoding
          else [decoded,m,s,in]= vitdec(depunctured,trel,Tbdepth,'cont','unquant',m,s,in);
         end
         decodeds = [decodeds  decoded]; 
      end % End of for i_sym  Loop
      range=delay+1:length(decodeds);  Lrange=length(range);
      nobe=nobe+sum(msgs(1:Lrange)~=decodeds(range)); notb=notb+Lrange;
      if nobe>Target_no_of_error,  break;  end
   end % End of for iter  Loop
   EbN0dB=EbN0dBs(i_EbN0);  BER= nobe/notb;  BERs(i_EbN0)=BER;  
   fprintf('\nWith EbN0dB=%3.1f, BER=%9.4e(=%d/%d)\n',EbN0dB,BER,nobe,notb)
end % End of for i_EbN0  Loop
BER_theory=prob_error(EbN0s_t+3,Mod,b,'BER');
subplot(222), semilogy(EbN0dBs,BERs,'r*', EbN0dBs_t,BER_theory,'b')