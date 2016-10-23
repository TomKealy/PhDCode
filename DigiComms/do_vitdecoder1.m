%do_vitdecoder1.m 
% shows various uses of Communication Toolbox function convenc()
% with KxN Code generator matrix Gc - octal polynomial representation 
clear, clf
%msg=[1 0 1 1 0 0 0];  
msg=randint(1,100)
lm=length(msg); % Message and its length
potbe=0.02; % Probability of transmitted bit error
Gc=[5 7]; % 1 0 1 -> 5, 1 1 1 -> 7 (octal number)
Lc=3; % 1xK constraint length vector for each input stream
[K,N]=size(Gc); % Number of encoder input/output bits
trel=poly2trellis(Lc,Gc); % Trellis structure
ch_input1=convenc(msg,trel); % Convolutional encoder
notbe1=ceil(potbe*length(ch_input1));
error_bits1=randerr(1,length(ch_input1),notbe1);
detected1= rem(ch_input1+error_bits1,2); % Received/modulated/detected 
% with hard decision
Tbdepth=3; % Traceback depth 
decoded1= vitdec(detected1,trel,Tbdepth,'trunc','hard')
noe_vitdec_trunc_hard=sum(msg~=decoded1(1:lm))
decoded2= vitdec(detected1,trel,Tbdepth,'cont','hard');
noe_vitdec_cont_hard=sum(msg(1:end-Tbdepth)~=decoded2(Tbdepth+1:end))
% with soft decision
ncode= [detected1+0.1*randn(1,length(detected1)) zeros(1,Tbdepth*N)];
quant_levels=[0.001,.1,.3,.5,.7,.9,.999];
NSDB=ceil(log2(length(quant_levels))); % Number of Soft Decision Bits
qcode= quantiz(ncode,quant_levels); % Quantized
decoded3= vitdec(qcode,trel,Tbdepth,'trunc','soft',NSDB);
noe_vitdec_trunc_soft=sum(msg~=decoded3(1:lm))
decoded4= vitdec(qcode,trel,Tbdepth,'cont','soft',NSDB);
noe_vitdec_cont_soft=sum(msg~=decoded4(Tbdepth+1:end))
% Repetitive use of vitdec() to process the data block by block
delay=Tbdepth*K; % Decoding delay depending on the traceback depth
% Initialize the message sequence, decoded sequence,
%   state metric, traceback state/input, and encoder state.
msg_seq=[]; decoded_seq=[]; 
m=[]; s=[]; in=[]; encoder_state=[];
N_Iter=100;
for itr=1:N_Iter
   msg=randint(1,1000); % Generate the message sequence in a random way
   msg_seq= [msg_seq msg]; % Accumulate the message sequence
   if itr==N_Iter, msg=[msg zeros(1,delay)]; end % Append with zeros
   [coded,encoder_state]=convenc(msg,trel,encoder_state);
   [decoded,m,s,in]=vitdec(coded,trel,Tbdepth,'cont','hard',m,s,in);
   decoded_seq=[decoded_seq decoded];
end
lm=length(msg_seq);
noe_repeated_use=sum(msg_seq(1:lm)~=decoded_seq(delay+[1:lm]))