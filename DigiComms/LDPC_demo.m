%LDPC_demo.m
clear 
M=5; N=10; K=N-M; % Length of message vector 
%H=gen_ldpc(M,N); % Generate a parity-check matrix
H=[ 1  1  0  1  1  0  0  1  1  1;
    1  0  1  1  1  0  1  1  1  0;
    0  1  1  1  0  1  1  1  0  1;
    0  1  0  0  1  1  1  0  0  1;
    1  0  1  0  0  1  0  0  1  0];
H1=H(:,1:N-K); H2=H(:,N-K+1:N);  H1T=H1.'; H2T=H2.';
Rc = K/N; a = 1;  % Code rate and Fading amplitude
EbN0dBs = [1 3 5];  N_EbN0dBs = length(EbN0dBs);
Nframes = 500; 
for nEN = 1:N_EbN0dBs
   EbN0=10^(EbN0dBs(nEN)/10); % convert Eb/N0[dB] to normal numbers
   EbN0r2=2*EbN0*Rc; EbN0r4=2*EbN0r2*a;  sigma=1/sqrt(EbN0r2);
   ner(nEN) = 0;
   for nframe = 1:Nframes
      m = randint(1,K,2);  % information message bits
      p = rem(m*rem(H2T*inv_GF2(H1T),2),2); % parity vector Eq.(9.4.50)
      c = [p m]; % code vector
      noise = sigma*randn(1,length(c));
      r = a*(2*c-1) + noise; % BPSK modulation and noise
      Lcy = EbN0r4.*r;
      y = LDPC_decoder(H,Lcy); % decoding
      ner(nEN) = ner(nEN) + sum(m~=y(end-K+1:end));
   end   %nframe
end   %nEN
ber = ner/(K*Nframes)
