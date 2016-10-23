function [output,state]=conv_encoder(G,K,input,state,termmode)
% generates the output sequence of a binary convolutional encoder
% G    : N x LK Generator matrix of a convolutional code
% K    : Number of input bits entering the encoder at each clock cycle.
% input: Binary input sequence
% state: State of the convolutional encoder
% termmode='trunc' for no termination with all-0 state 
%Copyleft: Won Y. Yang, wyyang53@hanmail.net, CAU for academic use only
if isempty(G), output=input; return; end
tmp= rem(length(input),K);
input=[input zeros(1,(K-tmp)*(tmp>0))];
[N,LK]=size(G);  
if rem(LK,K)>0
  error('The number of column of G must be a multiple of K!')
end
%L=LK/K;
if nargin<4|(nargin<5 & isnumeric(state))
  input= [input zeros(1,LK)]; %input= [input zeros(1,LK-K)]; end
end
if nargin<4|~isnumeric(state)
  state=zeros(1,LK-K); 
end
input_length= length(input);
N_msgsymbol= input_length/K;
input1= reshape(input,K,N_msgsymbol);
output=[]; 
for l=1:N_msgsymbol % Convolution output=G*input
   ub= input1(:,l).';
   [state,yb]= state_eq(state,ub,G);
   output= [output yb];
end