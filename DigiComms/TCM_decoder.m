function decoded_seq=TCM_decoder(state_eq,K,Nsb,received,Constellation,opmode)
% Performs the Viterbi algorithm on the PSK demodulated signal
% Input: state_eq = External function for state eqn saved in an M-file
%        K   = Number of input bits entering the encoder at each cycle.
%        Nsb = Number of state bits of the TCM encoder
%        received = received sequence
%        Constellation: Signal sets for signal mapper
%Copyleft: Won Y. Yang, wyyang53@hanmail.net, CAU for academic use only
N_states=2^Nsb;
N_msgsymbol=length(received);
for m=1:N_states
   for n=1:N_msgsymbol+1
      states(m,n)=0; %inactive in the trellis diagram
      p_state(m,n)=0; n_state(m,n)=0; input(m,n)=0;
   end   
end      
states(1,1)=1; % make the initial state active     
cost(1,1)=0; K2=2^K;
for n=1:N_msgsymbol
   y=received(n); %received sequence
   n1=n+1;
   for m=1:N_states
      if states(m,n)==1 %active
        xb=deci2bin1(m-1,Nsb);
        for m0=1:K2
           u=deci2bin1(m0-1,K);
           [nxb(m0,:),yb(m0)]=feval(state_eq,xb,u,Constellation);
           nxm0=bin2deci(nxb(m0,:))+1;  states(nxm0,n1)=1;
           % Accumulated squared Euclidean distance as path-to-node cost 
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
if nargin<6 | ~strncmp(opmode,'term',4)
  % trace back from the best-metric state (default)
  [min_cost,m]=min(cost(:,n1));
 else m=1; % trace back from the all-0 state
end    
for n=n1:-1:2
  decoded_seq= [deci2bin1(input(m,n),K) decoded_seq];
  m=p_state(m,n);
end