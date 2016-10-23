%Simple PLL m-file demonstration
%Run program from editor Debug (F5)
%JC 5/17/09
%This m-file demonstrates a PLL which tracks and demodulates an FM carrier.

clear all; 
close all; 

f=1000;%Carrier frequency 
fs=100000;%Sample frequency
N=5000;%Number of samples
Ts=1/fs;
t=(0:Ts:(N*Ts)- Ts);

%Create the message signal
f1=100;%Modulating frequency
msg=sin(2*pi*f1*t);
kf=.0628;%Modulation index

%Create the real and imaginary parts of a CW modulated carrier to be tracked.
Signal=exp(j*(2*pi*f*t+2*pi*kf*cumsum(msg)));%Modulated carrier
Signal1=exp(j*(2*pi*f*t));%Unmodulated carrier

%Initilize PLL Loop 
phi_hat(1)=30; 
e(1)=0; 
phd_output(1)=0; 
vco(1)=0; 

%Define Loop Filter parameters(Sets damping)
kp=0.15; %Proportional constant 
ki=0.1; %Integrator constant 
%PLL implementation 

for n=2:length(Signal) 
  vco(n)=conj(exp(j*(2*pi*n*f/fs+phi_hat(n-1))));%Compute VCO 
  
  phd_output(n)=imag(Signal(n)*vco(n));%Complex multiply VCO x Signal input 
  
  e(n)=e(n-1)+(kp+ki)*phd_output(n)-ki*phd_output(n-1);%Filter integrator 
  
  phi_hat(n)=phi_hat(n-1)+e(n);%Update VCO 
end; 

%Plot waveforms 
startplot = 1;
endplot = 1000;

figure(1);
subplot(3,2,1);
plot(t(startplot:endplot), msg(startplot:endplot));
title('100 Hz message signal');
%xlabel('Time (seconds)');
ylabel('Amplitude');
grid;

figure(1);
subplot(3,2,2);
plot(t(startplot:endplot), real(Signal(startplot:endplot)));
title('FM (1KHz carrier modulated with a 100 Hz message signal)');
%xlabel('Time (seconds)');
ylabel('Amplitude');
grid;

figure(1)
subplot(3,2,3);
plot(t(startplot:endplot), e(startplot:endplot));
title('PLL Loop Filter/Integrator Output');
%xlabel('Time (seconds)');
ylabel('Amplitude');
grid;

subplot(3,2,4);
plot(t(startplot:endplot), real(vco(startplot:endplot)));
title('VCO Output (PLL tracking the input signal)');
%xlabel('Time (seconds)');
ylabel('Amplitude');
grid;

subplot(3,2,5);
plot(t(startplot:endplot), phd_output(startplot:endplot));
title('Phase Detecter Output');
xlabel('Time (seconds)');
ylabel('Amplitude');
grid;

subplot(3,2,6);
plot(t(startplot:endplot), real(Signal1(startplot:endplot)));
title('Unmodulated Carrier');
xlabel('Time (seconds)');
ylabel('Amplitude');
grid;