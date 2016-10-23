%do_sync_for_DMB.m
clear, clf
Nfft=2048; Ng=504; Nnull=2656;
CFO = -1.7;  phase = 0; % A pseudo Carrier Frequency/Phase Offset
nF = 1 % A pseudo Frame Time Offset (delay)
PRS = phase_ref_symbol;
prs = ifft(PRS,Nfft); % Time-domain PRS 
tx = [prs(Nfft-Ng+1:Nfft) prs]; % Add Cyclic Prefix
L_tx=length(tx);
% Set up the (pseudo) CFO (pretended)
Nd = 100; % Delay tolerance
y = set_CFO(tx,CFO,phase,Nfft);  A=0.05; 
y = [A*(rand(1,Nd)-0.5+j*(rand(1,Nd)-0.5)) y];
y = [y A*(rand(1,Nd)-0.5+j*(rand(1,Nd)-0.5))];
if nF>0
  y=[A*(rand(1,nF)-0.5+j*(rand(1,nF)-0.5)) y(1:end-nF)]; %Delayed
 elseif nF<0
  y=[y(1-nF:end) A*(rand(1,-nF)-0.5+j*(rand(1,-nF)-0.5))]; %Advanced
end     
% IFO (Integral Frequency Offset) estimation
IFO_range = [-3:3];
y_ifft = y(Nd+Ng+[1:Nfft]);
IFO_est = IFO_estimate(y_ifft,PRS,IFO_range); % Eq.(11.3.9a,b,c)
%In order to realize how critical the IFO estimate is for FFS, 
% activate the following statement even with a very short delay nF=1;
% IFO_est = zeros(1,3);
% FFS (fine frame synchronization)
Y = fft(y_ifft,Nfft);
YX = Y.*conj(rotate_r(PRS,IFO_est));
[Max,nF_] = max(abs(ifft(YX)));  % Eq.(11.3.11)
if nF_>Nfft/2,  nF_=nF_-Nfft;  end % Periodicity of FFT/IFFT
nF_h = (nF_-1)
% FFO (Fractional Frequency Offset) estimation
nn = Nd+[1:Ng]; nn1 = Nfft + nn;
% To realize the importance of reflecting the FFS into FFO estimation, 
% activate the following statement with a delay nF=2;
% nF_h = 0;
FFO_est = angle(y(nF_h+nn1)*y(nF_h+nn)')/(2*pi) % Eq.(11.3.12)
% FFO estimate without incorporating the FFS result
FFO_est0 = angle(y(nn1)*y(nn)')/(2*pi) % Eq.(11.3.12) without nF_h
% Overall CFO estimate
CFO_est = IFO_est + FFO_est;
fprintf('\n For CFO=%10.8f,\n', CFO);
form = '  CFO_estimate(%10.8f) = IFO(%10.8f)+FFO(%10.8f)\n';
for i=1:3
   fprintf(form, CFO_est(i), IFO_est(i), FFO_est);
end   
% CFO compensated symbols
y_compensated = compensate_CFO(y,CFO_est(3),Nfft);
% Discrepancy between the CFO compensated symbols and original ones
discrepancy = norm(tx-y_compensated(Nd+[1:L_tx]))/L_tx