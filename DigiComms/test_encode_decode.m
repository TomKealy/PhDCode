%test_encode_decode.m  to try using encode()/decode()
N=7; K=4; % Codeword (Block) length and Message size
g=cyclpoly(N,K); % Generator polynomial for a cyclic (N,K) code
Nm=10; % # of K-bit message vectors
msg=randint(Nm,K); % Nm x K message matrix 
coded = encode(msg,N,K,'cyclic',g); % Encoding
% Add bit errors with transmitted BER potbe=0.1
potbe=0.1; received=rem(coded+randerr(Nm,N,[0 1;1-potbe potbe]),2); 
decoded=decode(received,N,K,'cyclic',g); % Decoding
% Probability of message bit errors after decoding/correction
pobe=sum(sum(decoded~=msg))/(Nm*K) % BER
% Usage of rsenc()/rsdec()
M=3; % Galois Field integer corresponding to the # of bits per symbol
N=2^M-1; K=3; dc=(N-K)/2; % Codeword length and Message size 
msg=gf(randint(Nm,K,2^M),M); % Nm x K GF(2^M) Galois Field msg matrix  
coded = rsenc(msg,N,K); % Encoding
noise = randerr(Nm,N,[1 dc+1]).*randint(Nm,N,2^M);
received = coded+noise; % Add a noise
[decoded,numerr]=rsdec(received,N,K); % Decoding
[msg decoded], numerr, pose=sum(sum(decoded~=msg))/(Nm*K) % SER
