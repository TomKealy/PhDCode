%do_vitdecoder.m 
% Try using conv_encoder()/vit_decoder()
clear, clf
msg=[1 0 1 1 0 0 0];  % msg=randint(1,100)
lm=length(msg); % Message and its length
G=[1 0 1;1 1 1]; % N x LK Generator polynomial matrix
K=1; N=size(G,1); % Size of encoder input/output
potbe=0.02; % Probability of transmitted bit error
% Use of conv_encoder()/vit_decoder()
ch_input=conv_encoder(G,K,msg) % Self-made convolutional encoder
notbe=ceil(potbe*length(ch_input));
error_bits=randerr(1,length(ch_input),notbe); 
detected= rem(ch_input+error_bits,2); % Received/modulated/detected
decoded= vit_decoder(G,K,detected)
noe_vit_decoder=sum(msg~=decoded(1:lm))  % Number of errors