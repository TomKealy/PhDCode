%dc09p07.m
% To practice using convenc() and vitdec() for channel coding
clear, clf
Gc=[4 5 11;1 4 2]; % Octal code generator matrix
K=size(Gc,1); % Number of encoder input bits
% Constraint length vector 
Gc_m=max(Gc.');
for i=1:length(Gc_m), Lc(i)=length(deci2bin1(oct2dec(Gc_m(i))));  end
trel=poly2trellis(Lc,Gc);
Tbdepth=sum(Lc)*5; delay=Tbdepth*K;
lm=1e5; msg=randint(1,lm);
transmission_ber=0.02;
notbe=round(transmission_ber*lm); % Number of transmitted bit errors
ch_input=convenc([msg zeros(1,delay)],trel); 
% Received/modulated/detected signal
ch_output= rem(randerr(1,length(ch_input),notbe)+ch_input,2);
decoded_trunc= vitdec(ch_output,trel,Tbdepth,'trunc','hard');
ber_trunc= sum(msg~=decoded_trunc(????))/lm;
decoded_cont= vitdec(ch_output,trel,Tbdepth,'cont','hard');
ber_cont=sum(msg~=decoded_cont(????????????))/lm;
% It is indispensable to use the delay for the decoding result 
%  obtained using vitdec(,,,'cont',) 
nn=[0:100-1]; 
subplot(221), stem(nn,msg(nn+1)), title('Message sequence')
subplot(223), stem(nn,decoded_cont(nn+1)), hold on 
stem(delay,0,'rx')
decoded_term= vitdec(ch_output,trel,Tbdepth,'term','hard');
ber_term=sum(msg~=decoded_term(????))/lm;
fprintf('\n BER_trunc  BER_cont   BER_term')
fprintf('\n %9.2e  %9.2e  %9.2e\n', ber_trunc,ber_cont,ber_term)