%do_cyclic_codes.m
clear, clf
gsymbols=['bo';'r+';'kx';'md';'g*'];
% Theoretical BER curves
SNRdBs=0:0.01:14;  SNRs=10.^(SNRdBs/10); pemb_uncoded=Q(sqrt(SNRs)); 
semilogy(SNRdBs,pemb_uncoded,':'), hold on
use_decode = 0; % Use encode()/decode or not
MaxIter=1e5; Target_no_of_error=100;  SNRdBs1=2:4:10;
for iter=1:2
   if iter==1, N=15; K=7; nceb=2; %(15,7) BCH code generator 
    else    N=31; K=16; nceb=3; % (31,16) BCH code generator
   end
   gBCH=bchgenpoly(N,K); %Galois row vector representing (N,K) BCH code
   g=double(gBCH.x); 
   clear('S')
   if use_decode>0, H=cyclgen(N,g); E=syndtable(?);
    else
     E= combis(N,nceb); % All error patterns
     for i=1:size(E,1)
        S(i,:) = cyclic_decoder0(E(i,:),N,K,g); % Syndrome
        epi(bin2deci(S(i,:))) = ?; % Error pattern indices
     end
   end
   Rc = K/N;  %code rate
   SNRcs=SNRs*Rc; sqrtSNRcs=sqrt(SNRcs);
   ets=Q(sqrt(SNRcs)); % transmitted bit error probability Eq.(7.3.4)
   pemb_t=prob_err_msg_bit(ets,N,nceb);
   semilogy(SNRdBs,pemb_t,gsymbols(iter,1),'Markersize',5)
   for iter1=1:length(SNRdBs1)
      SNRdB = SNRdBs1(iter1); SNR = 10.^(SNRdB/10);   
      SNRc = SNR*Rc; sqrtSNRc = sqrt(SNRc); 
      et=Q(sqrt(SNRc)); % transmitted bit error probability Eq.(7.3.4)
      K100 = K*100;
      nombe = 0; 
      for iter2=1:MaxIter
         msg=randint(1,K100); % Message vector
         if use_decode==0,  coded=cyclic_???????(msg,N,K,g); 
          else  msg=reshape(msg,length(msg)/K,K);
                coded=??????(msg,N,K,'cyclic',g);
         end
         r= 2*coded-1 + randn(size(coded))/sqrtSNRc; % Received vector
         r_sliced= 1*(r>0); % Sliced
         if use_decode==0 
           decoded= cyclic_???????(r_sliced,N,K,g,E,epi);
          else  decoded= ??????(r_sliced,N,K,'cyclic',g,E);  
         end
         nombe= nombe + sum(sum(decoded~=msg)); 
         if nombe>Target_no_of_error, break; end
      end
      lm=iter2*K100; pemb=nombe/lm; % Message bit error probability
      semilogy(SNRdB,pemb,gsymbols(iter,:),'Markersize',5)
   end   
end