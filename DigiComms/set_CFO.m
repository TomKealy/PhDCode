function tx_with_CFO = set_CFO(tx,CFO,phase,Nfft,STO)
% CFO/phase = pseudo Carrier frequency/phase offset (pretended)
% STO       = Symbol Time Offset
if nargin<5, STO = 0; end % Symbol Time Offset
if nargin<4, Nfft = 64; end
if nargin<3, phase = 0; end
m = -STO + [0:length(tx)-1];
tx_with_CFO = tx.*exp(j*(2*pi*CFO/Nfft*m+phase));