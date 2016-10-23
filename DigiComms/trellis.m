function [nout,nstate,pout,pstate] = trellis(G)
% copyright 1998, Yufei Wu, MPRG lab, Virginia Tech for academic use
% set up the trellis with code generator G in binary matrix form. 
% G: Generator matrix with feedback/feedforward connection in row 1/2
% e.g. G=[1 1 1; 1 0 1] for the turbo encoder in Fig. 9.15(a) 
% nout(i,1:2): Next output [xs=m xp](-1/+1) for state=i, message in=0 
% nout(i,3:4): next output [xs=m xp](-1/+1) for state=i, message in=1
% nstate(i,1): next state(1,...2^M) for state=i, message input=0
% nstate(i,2): next state(1,...2^M) for state=i, message input=1
% pout(i,1:2): previous out [xs=m xp](-1/+1) for state=i, message in=0
% pout(i,3:4): previous out [xs=m xp](-1/+1) for state=i, message in=1
% pstate(i,1): previous state having come to state i with message in=0
% pstate(i,2): previous state having come to state i with message in=1
% See Fig. 9.16 for the meanings of the output arguments.
[N,L] = size(G); % Number of output bits and Consraint length
M=L-1;  Ns=2^M; % Number of bits per state and Number of states
% Set up next_out and next_state matrices for RSC code generator G
for state_i=1:Ns
   state_b = deci2bin1(state_i-1,M); % Binary state
   for input_bit=0:1
      d_k = input_bit; 
      a_k=rem(G(1,:)*[d_k  state_b]',2); % Feedback in Fig.9.15(a) out(input_bit+1,:)=[d_k  rem(G(2,:)*[a_k  state_b]',2)]; % Forward
      state(input_bit+1,:)=[a_k  state_b(1:M-1)]; % Shift register
   end
   nout(state_i,:) = 2*[out(1,:) out(2,:)]-1; % bipolarize
   nstate(state_i,:) = [bin2deci(state(1,:)) bin2deci(state(2,:))]+1;
end
% Possible previous states having reached the present state 
%  with input_bit=0/1
for input_bit=0:1
   bN = input_bit*N; b1 = input_bit+1; % Number of output bits = 2;
   for state_i=1:Ns
      pstate(nstate(state_i,b1),b1) = state_i;
      pout(nstate(state_i,b1),bN+[1:N]) = nout(state_i,bN+[1:N]);
   end 
end