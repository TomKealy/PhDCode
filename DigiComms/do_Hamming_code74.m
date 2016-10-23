function do_Hamming_code74(SNRbdB,MaxIter)
% (7,4) Hamming code
if nargin<2, MaxIter=1e6; end
n=3; N=2^n-1; K=2^n-1-n; % Codeword (Block) length and Message size
Rc=K/N; % Code rate
SNRb=10.^(SNRbdB/10); 
SNRbc=SNRb*Rc; 
sqrtSNRbc=sqrt(SNRbc);
pemb_uncoded=Q(sqrt(SNRb)); % Uncoded msg BER with BPSK by Eq.(7.3.4)
et=Q(sqrt(SNRbc)); % Crossover probability by Eq.(7.3.4)
L=2^K;
for i=1:L
   M(i,:)=deci2bin1(i-1,K); % All message vectors
end 
[H,G]=Hammgen(n); % [H,G]=Ham_gen(n): Eq.(9.4.12)&(9.4.14)
Hamming_code =rem(M*G,2); % Eq.(9.4.13)
Min_distance =min(sum(Hamming_code(2:L,:)')); % Eq.(9.4.8)
No_of_correctable_error_bits=floor((Min_distance-1)/2); % Eq.(9.4.10) 
E= combis(N,No_of_correctable_error_bits); % Error patterns (Fig.9.8)
S= rem(E*H',2); NS=size(S,1); % The syndrome matrix
nombe=0; 
for iter=1:MaxIter
   msg=randint(1,K); % Message vector
   coded= rem(msg*G,2); % Coded vector
   modulated= 2*coded-1; % BPSK-modulated vector
   r= modulated +randn(1,N)/sqrtSNRbc; % Received vector with noise
   r_sliced= r>0; % Sliced
   r_c=r_sliced; % To be corrected only if it has a syndrome
   s= rem(r_sliced*H',2); % Syndrome
   for m=1:NS  % Error correction depending on the syndrome
      if s==S(m,:), r_c=rem(r_sliced+E(m,:),2);  break;  end
   end
   nombe=nombe+sum(msg~=r_c(N-K+1:N));  if nombe>100, break; end
end
pemb=nombe/(K*iter); % Message bit error probability
pemb_t=prob_err_msg_bit(et,N,No_of_correctable_error_bits); 
fprintf('\n Uncoded Messsage BER=%8.6f',pemb_uncoded)
fprintf('\nMessage BER=%5.4f(Theoretical BER=%5.4f)\n',pemb,pemb_t)