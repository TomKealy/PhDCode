function CFO_compensated = compensate_CFO(tx,CFO_est,Nfft,STO)
% CFO_est = Estimated carrier frequency offset 
% STO     = Symbol Time Offset 
if nargin<4, STO = 0; end
if nargin<3, Nfft = 64; end
m = -STO + [0:length(tx)-1]; %m = -STO+m0 +[0:length(tx)-1];
CFO_compensated = tx.*exp(-j*2*pi*CFO_est/Nfft*m);