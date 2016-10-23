function coded= cyclic_encoder(msg_seq,N,K,g)
% Cyclic (N,K) encoding of input msg_seq m with generator polynomial g 
Lmsg=length(msg_seq); Nmsg=ceil(Lmsg/K);
Msg= [msg_seq(:); zeros(Nmsg*K-Lmsg,1)];  
Msg= reshape(Msg,K,Nmsg).';
coded= []; 
for n=1:Nmsg
   msg= Msg(n,:);
   for i=1:N-K, x(i)=0; end
   for k=1:K
      tmp= rem(msg(K+1-k)+x(N-K),2); % msg(K+1-k)+g(N-K+1)*x(N-K)
      for i=N-K:-1:2,  x(i)= rem(x(i-1)+g(i)*tmp,2);  end
      x(1)=g(1)*tmp;
   end
   coded= [coded x msg]; % Eq.(9.4.26)
end