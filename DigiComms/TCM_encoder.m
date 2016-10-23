function [output,state]=TCM_encoder(state_eq,K,Ns,N,input,state,Constellation,opmode)
% Generates the output sequence of a binary TCM encoder
% Input: state_eq = External function for state eqn saved in an M-file
%        K    = Number of input bits entering the encoder at each cycle
%        Nsb  = Number of state bits of the TCM encoder
%        N    = Number of output bits of the TCM encoder
%        input= Binary input message seq.
%        state= State of the conv_encoder
%        Constellation= Signal sets for signal mapper
% Output: output= Sequence of signal points on Contellation
%         state = Updated state
%Copyleft: Won Y. Yang, wyyang53@hanmail.net, CAU for academic use only
tmp= rem(length(input),K); 
input= [input zeros(1,(K-tmp)*(tmp>0))];
if nargin<6, state=zeros(1,Nsb); 
 elseif length(state)==2^N, Constellation=state; state=zeros(1,Nsb);
end
input_length= length(input);  N_msgsymbol= input_length/K;
input1= reshape(input,K,N_msgsymbol).';
outputs= []; 
for l=1:N_msgsymbol
   ub= input1(l,:);
   [state,output]=feval(state_eq,state,ub,Constellation);
   outputs= [outputs output];
end
