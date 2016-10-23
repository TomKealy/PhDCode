%do_CFO.m
clear, clf
Nfft = 64; Ng = 16;
CFO = 1.7; phase = 0; STO = -1; % CFO/Phase Offset/STO 
NB = 6; % Which block to start computing the correlation with? 
Nw = Ng; % Correlation window size
[short_preamble,S] = short_train_seq(Nfft);
[long_preamble,L] = long_train_seq(Nfft);
% Time-domain training symbol
tx = [short_preamble  long_preamble]; 
% Set up a pseudo CFO
tx_offset = set_CFO(tx,CFO,phase,Nfft); 
% Coarse CFO estimation
coarse_CFO = coarse_CFO_estimate(tx_offset(1:160),NB,Nw,STO);
% Fine CFO estimation
fine_CFO = fine_CFO_estimate(tx_offset(161:320),coarse_CFO,Nfft);
% Overall CFO estimate
CFO_estimate = coarse_CFO + fine_CFO;
form1 = '\n For CFO=%10.8f, CFO estimate(%10.8f) = ';
form2 = 'coarse estimate(%10.8f)+fine estimate(%10.8f)\n';
fprintf([form1 form2],CFO,CFO_estimate,coarse_CFO,fine_CFO);
% CFO compensated symbols
tx_compensated = compensate_CFO(tx_offset,CFO_estimate,Nfft,STO);
discrepancy = norm(tx-tx_compensated)/length(tx)