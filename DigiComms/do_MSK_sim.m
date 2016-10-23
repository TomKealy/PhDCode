%do_MSK_sim.m
clear, clf
SNRdBt=[0:0.1:10]; SNRdBs=[5  10]; EbN0dBs= SNRdBs-3;
b=1; M=2^b; % Number of bits per symbol and Modulation order
Tb=1; Nb=16; T=Tb/Nb; Ts=b*Tb; % Bit duration and Sample time
wc=4*pi/Tb; % Carrier frequency[rad/s]
t=[0:2*Nb-1]*T;  pihTbt=pi/2/Tb*t;
suT=sqrt(2/Tb)*T*[cos(pihTbt-pi/2).*cos(wc*(t-Tb)); sin(pihTbt).*sin(wc*t)];  % Eq. (7.6.8a,b)
M_filter1=fliplr(suT(1,:)); M_filter2=fliplr(suT(2,:));
Target_no_of_error = 50;
Simulink_mdl='MSK_passband_sim'; 
for i=1:length(SNRdBs)
   SNRdB = SNRdBs(i);  EbN0dB = SNRdB-3;
   sim(Simulink_mdl,1e5*Tb),  BERs(i) = ber(end,1);
end
pobet_PSK=prob_error(SNRdBt,'PSK',b,'bit'); 
pobet_FSK=prob_error(SNRdBt,'FSK',b,'bit');
semilogy(SNRdBt,pobet_PSK,'k', SNRdBt,pobet_FSK,'k:', SNRdBs,BERs,'b*')
