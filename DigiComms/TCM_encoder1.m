function [outputs,state]=TCM_encoder1(Gc,input,state,Constellation)
% generates the output sequence of a binary TCM encoder
% Input:  Gc    = Code generator matrix consisting of octal numbers
%         K     = Number of input bits entering the encoder 
%         Nsb   = Number of state bits of the TCM encoder
%         N     = Number of output bits of the TCM encoder
%         input = Binary input message sequence
%         state = State of the TCM encoder
%         Constellation = Signal sets for signal mapper
% Output: output = a sequence of signal points on Contellation
% Copyleft: Won Y. Yang, wyyang53@hanmail.net, CAU for academic use only
[K,N]=size(Gc);
for k=1:K, nsb(k)=length(de2bi(oct2dec(max(Gc(k,:)))))-1; end
Nsb=sum(nsb);   
tmp= rem(length(input),K); 
input= [input zeros(1,(K-tmp)*(tmp>0))];
if nargin<3, state=zeros(1,Nsb); 
  elseif length(state)==2^N
    Constellation=state; state=zeros(1,Nsb);
end
input_length= length(input);
N_msgsymbol= input_length/K;
input1= reshape(input,K,N_msgsymbol).';
outputs= [];
for l=1:N_msgsymbol
   ub= input1(l,:);
   is=1; nstate=[]; output=zeros(1,N);
   for k=1:K
      tmp = [ub(k) state(is:is+nsb(k)-1)]; 
      nstate = [nstate tmp(1:nsb(k))]; is=is+nsb(k); 
      for i=1:N
         output(i) = output(i) + deci2bin1(oct2dec(Gc(k,i)),nsb(k)+1)*tmp';
      end
   end
   state = nstate; 
   output = Constellation(bin2deci(rem(output,2))+1);
   outputs = [outputs output];
end
