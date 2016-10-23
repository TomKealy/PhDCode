function decoded_seq=TCM_decoder1(Gc,demod,Constellation,opmode)
% performs the Viterbi algorithm on the PSK demodulated signal
% Inpit: state_eq = External function for state equation saved
%                    in an M-file
%        K   = Number of input bits entering encoder at each cycle.
%        Nsb = Number of state bits of TCM encoder
%Copyleft: Won Y. Yang, wyyang53@hanmail.net, CAU for academic use only
[K,N]=size(Gc);
for k=1:K, nsb(k)=length(de2bi(oct2dec(max(Gc(k,:)))))-1; end
Nsb=sum(nsb); 
N_states=2^Nsb;
N_msgsymbol=length(demod);
for m=1:N_states
   for n=1:N_msgsymbol+1
      states(m,n)=0; % inactive in the trellis diagram
      p_state(m,n)=0; n_state(m,n)=0; input(m,n)=0;
   end   
end      
states(1,1)=1; % make the initial state active     
cost(1,1)=0; K2=2^K;
for n=1:N_msgsymbol
   y=demod(n); %received sequence
   n1=n+1;
   for m=1:N_states
      if states(m,n)==1 %active
        xb=deci2bin1(m-1,Nsb);
        for m0=1:K2
           u=deci2bin1(m0-1,K);
           is=1; nstate=[]; output=zeros(1,N);
           for k=1:K
              tmp = [u(k) xb(is:is+nsb(k)-1)]; 
              nstate = [nstate tmp(1:nsb(k))]; is=is+nsb(k); 
              for i=1:N
                 output(i) = output(i) + deci2bin1(oct2dec(Gc(k,i)),nsb(k)+1)*tmp';
              end
           end
           nxb(m0,:) = nstate; 
           yb(m0) = Constellation(bin2deci(rem(output,2))+1);
           nxm0=bin2deci(nxb(m0,:))+1;
           states(nxm0,n1)=1;
           difference=y-yb(m0);
           d(m0)=cost(m,n)+difference*conj(difference);
           if p_state(nxm0,n1)==0 
             cost(nxm0,n1)=d(m0);
             p_state(nxm0,n1)=m; input(nxm0,n1)=m0-1; 
            else
             [cost(nxm0,n1),i]=min([d(m0) cost(nxm0,n1)]);
             if i==1, p_state(nxm0,n1)=m; input(nxm0,n1)=m0-1; end
           end 
        end        
      end
   end
end
decoded_seq=[];
if nargin<4 | ~strncmp(opmode,'term',4)
  % trace back from the best-metric state (default)
  [min_cost,m]=min(cost(:,n1)); 
 else m=1; %trace back from the all-0 state
end    
for n=n1:-1:2
   decoded_seq= [deci2bin1(input(m,n),K) decoded_seq];
   m=p_state(m,n);
end
