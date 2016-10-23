%do_rcos2.m
clear, clf
r=0.5; % Roll-off factor
Ts=1; T=Ts/4; % Symbol time and Sample time of filter input/output
Fd=1/Ts; Fs=Fd*(Ts/T); %Sampling frequencies of filter input/output
Delay=1; N_data=20; x = randsrc(1,N_data);  
t_x = [0:N_data-1]*Ts; t_xd = t_x+Delay*Ts;
subplot(311), stem(t_x,x,'k'), axis([0 20 -1.3 1.3])
% Design RC filter.
r=0.5; RC_filter1 = rcosine(Fd,Fs,'fir',r,Delay);
% Filtering with upsampling
[yRC,t_yRC] = rcosflt(x,Fd,Fs,'filter',RC_filter1);
% One-shot method to design and filter data 
%  by using rcosflt() without "filter" in the 4th input argument 
[yRC1,t_yRC1] = rcosflt(x,Fs,Fs,'fir',r,Delay);
subplot(312), stem(t_xd,x) % Input signal delayed by Delay
hold on, plot(t_yRC,yRC,'bo', t_yRC1,yRC1,'bx')
axis([0 20 -1.3 1.3])
% Design a square-root raised-cosine (SRRC) filter.
SRRC_filter = rcosine(Fd,Fs,'fir/sqrt',r,Delay);
% SRRC Filtering with upsampling at XMTR
[ySRRC1,t_ySRRC1] = rcosflt(x,Fd,Fs,'filter',SRRC_filter);
% the same SRRC Filtering without upsampling at RCVR
[ySRRC2,t_ySRRC2] = rcosflt(ySRRC1,Fs,Fs,'filter/Fs',SRRC_filter); 
% where "/Fs" is used to filter without upsampling.
subplot(313), stem(t_xd,x) % the input signal delayed by Delay
hold on, plot(t_ySRRC1,ySRRC1,'.', t_ySRRC2,ySRRC2,'ro')
axis([0 20 -1.3 1.3])
