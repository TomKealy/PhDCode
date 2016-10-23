%do_cyclic_code.m 
% tries with a cyclic code.
clear
N=7; K=4; % N=15; K=7; % Codeword (Block) length and Message size
%N=31; K=16;
g=cyclpoly(N,K); g_=fliplr(g);
lm=5*K; msg= randint(1,lm);
% N=7; K=4; g_=[1 1 0 1]; g=fliplr(g_); msg=[1 0 1 1];
coded = cyclic_encoder(msg,N,K,g_); lc=length(coded);
no_transmitted_bit_errors=ceil(lc*0.05);
errors=randerr(1,lc,no_transmitted_bit_errors);
r = rem(coded+errors,2); % Received sequence
decoded = cyclic_decoder(r,N,K,g_); nobe=sum(decoded~=msg) 
coded1 = encode(msg,N,K,'cyclic',g); % Use the Communication Toolbox 
r1 = rem(coded1+errors,2); % Received sequence
decoded1 = decode(r1,N,K,'cyclic',g); nobe1=sum(decoded1~=msg)