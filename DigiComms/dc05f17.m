%dc05f17.m 
% plots the bit error probability vs. SNR in Figs. 5.13 and 5.17
clear, clf
SNRdBs=0:0.2:20; SNRbs=10.^(SNRdBs/10); %SNRt=Eav,b/(N0/2)
sqpi=sqrt(pi); tol=1e-6; % tolerance used for 'quad' integration
for b=1:4
   M=2^b;  tmp1=2*(M-1)/M/b; tmp2=(M/2)/(M-1); 
   poe_m_level= tmp1*Q(sqrt(3*b*SNRbs/(M*M-1))); %Eq.(5.2.43)
   for m=1:length(SNRbs) % Eq.(5.2.52)
      SNRb=SNRbs(m); 
      poe_m_dim(m)= tmp2*(1-Gauss_Hermite('f5252',10,SNRb,b)/sqpi); 
      %poe_m_dim(m)=tmp2*(1-quad('f5252_0',-20,20,tol,[],SNRb,b)/sqpi);
   end  
   semilogy(SNRdBs,poe_m_dim,'b',SNRdBs,poe_m_level,'r:'), hold on
end
