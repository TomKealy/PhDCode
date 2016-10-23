function x=cyclic_decoder0(r,N,K,g)
% Cyclic (N,K) decoding of an N-bit code r with generator polynomial g
for i=1:N-K, x(i)=r(i+K); end
for n=1:K
   tmp=x(N-K);
   for i=N-K:-1:2,  x(i)=rem(x(i-1)+g(i)*tmp,2);   end
   x(1)=rem(g(1)*tmp+r(K+1-n),2);
end