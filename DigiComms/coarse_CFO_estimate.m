function coarse_CFO_est = coarse_CFO_estimate(tx,NB,Nw,STO)
% Input:
%    tx = Received signal
%    NB = Which block of size Nw to start computing correlation with?
%    Nw = Correlation window size
%    STO= Symbol Time Offset 
% Output: coarse_CFO_est = Estimated carrier frequency offset 
if nargin<4, STO = 0; end
if nargin<3, Nw = 16; end
if nargin<2, NB = 6; end
for i=1:2
   nn = STO + Nw*(NB+i) + [1:Nw];
   CFO_est(i) = angle(tx(nn+Nw)*tx(nn)'); % Eq.(11.3.2a)
end   
CFO_est = sum(CFO_est)/pi; % Average
coarse_CFO_est = CFO_est - mod(CFO_est,4/128); % Stored with 8 bits