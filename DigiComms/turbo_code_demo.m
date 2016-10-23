%turbo_code_demo.m
% simulates the classical turbo encoding-decoding system.
% 1st encoder is terminated with tails bits. (lm+M) bits are scrambled 
% and passed to 2nd encoder, which is left open without termination.
clear
dec_alg = 1; % 0/1 for Log-MAP/SOVA
puncture = 1; % puncture or not 
rate = 1/(3-puncture);  % Code rate 
lu = 1000; % Frame size
Nframes = 100; % Number of frames
Niter = 4; % Number of iterations
EbN0dBs = 2.6; %[1 2 3];
N_EbN0dBs = length(EbN0dBs);
G = [1 1 1; 1 0 1]; % Code generator
a = 1; % Fading amplitude; a=1 in AWGN channel
[N,L]=size(G);  M=L-1;  lm=lu-M; % Length of message bit sequence
for nENDB = 1:N_EbN0dBs
   EbN0 = 10^(EbN0dBs(nENDB)/10); % convert Eb/N0[dB] to normal number
   L_c = 4*a*EbN0*rate;         % reliability value of the channel
   sigma = 1/sqrt(2*rate*EbN0); % standard deviation of AWGN noise
   noes(nENDB,:) = zeros(1,Niter);
   for nframe = 1:Nframes
      m = round(rand(1,lm));   % information message bits
      [temp,map] = sort(rand(1,lu));  % random interleaver mapping
      x = encoderm(m,G,map,puncture); % encoder output [x(+1/-1)
      noise = sigma*randn(1,lu*(3-puncture));
      r = a.*x + noise; % received bits
      y = demultiplex(r,map,puncture); % input for decoder 1 and 2
      Ly = 0.5*L_c*y; % Scale the received bits
      for iter = 1:Niter
         % Decoder 1
         if iter<2, Lu1=zeros(1,lu); % Initialize extrinsic information      
          else Lu1(map)=L_e2; % (deinterleaved) a priori information 
         end
         if dec_alg==0, L_A1=logmap(Ly(1,:),G,Lu1,1); % all information
          else          L_A1=sova(Ly(1,:),G,Lu1,1); % all information
         end
         L_e1= L_A1-2*Ly(1,1:2:2*lu)-Lu1; % Eq.(9.4.47)
         % Decoder 2         
         Lu2 = L_e1(map); % (interleaved) a priori information
         if dec_alg==0, L_A2=logmap(Ly(2,:),G,Lu2,2); % all information  
          else          L_A2=sova(Ly(2,:),G,Lu2,2); % all information 
         end
         L_e2= L_A2-2*Ly(2,1:2:2*lu)-Lu2; % Eq.(9.4.47)
         mhat(map)=(sign(L_A2)+1)/2; % Estimate the message bits
         noe(iter)=sum(mhat(1:lu-M)~=m); % Number of bit errors
      end  % End of iter loop
      % Total number of bit errors for all iterations
      noes(nENDB,:) = noes(nENDB,:) + noe;
      ber(nENDB,:) = noes(nENDB,:)/nframe/(lu-M); % Bit error rate
      fprintf('\n');
      for i=1:Niter, fprintf('%14.4e ', ber(nENDB,i)); end
   end   % End of nframe loop
end   % End of nENDB loop