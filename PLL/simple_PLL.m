%Simple PLL m-file demonstration
%Run program from editor Debug (F5)
%JC 5/17/09
%This m-file demonstrates a PLL which tracks and demodulates an FM carrier.
clear all; 
close all; 

fc = 1000; %Carrier frequency 
fs = 100000; %Sample frequency
N = 5000; %Number of samples
Ts = 1/fs;
t = (0:Ts:(N*Ts)- Ts);

%Create the message signal
f1 = 100; %Modulating frequency
msg = sin(2*pi*f1*t);
kf = .0628; %Modulation index

%Create the real and imaginary parts of a CW modulated carrier to be tracked.
Signal_carrier = exp(1j*(2*pi*fc*t+2*pi*kf*cumsum(msg))); %Modulated carrier
Signal1 = exp(1j*(2*pi*fc*t)); %Unmodulated carrier

%Initilize PLL Loop 
phi_hat(1) = 30; 
e(1) = 0; 
phase_dectector_output(1) = 0; 
vco(1) = 0; 

%Define Loop Filter parameters(Sets damping)
proportional_constant = 0.15; %Proportional constant 
integrator_constant = 0.1; %Integrator constant 

%PLL implementation 
for n=2:length(Signal_carrier) 
    vco(n)=conj(exp(1j*(2*pi*n*fc/fs+phi_hat(n-1)))); %Compute VCO 
    phase_dectector_output(n)=imag(Signal_carrier(n)*vco(n)); %Complex multiply VCO x Signal input 
    e(n)=e(n-1)+(proportional_constant+integrator_constant)*phase_dectector_output(n)-integrator_constant*phase_dectector_output(n-1); %Filter integrator 
    phi_hat(n)=phi_hat(n-1)+e(n); %Update VCO 
end;

%Plot waveforms 
startplot = 1;
endplot = 1000;

figure(1);
subplot(3, 2, 1);
plot(t(startplot:endplot), msg(startplot:endplot));
title('100 Hz message signal');
%xlabel('Time (seconds)');
ylabel('Amplitude');
grid;

figure(1);
subplot(3, 2, 2);
plot(t(startplot:endplot), real(Signal_carrier(startplot:endplot)));
title('FM (1KHz carrier modulated with a 100 Hz message signal)');
%xlabel('Time (seconds)');
ylabel('Amplitude');
grid;

figure(1)
subplot(3, 2, 3);
plot(t(startplot:endplot), e(startplot:endplot));
title('PLL Loop Filter/Integrator Output');
%xlabel('Time (seconds)');
ylabel('Amplitude');
grid;

subplot(3, 2, 4);
plot(t(startplot:endplot), real(vco(startplot:endplot)));
title('VCO Output (PLL tracking the input signal)');
%xlabel('Time (seconds)');
ylabel('Amplitude');
grid;

subplot(3, 2, 5);
plot(t(startplot:endplot), phase_dectector_output(startplot:endplot));
title('Phase Detecter Output');
xlabel('Time (seconds)');
ylabel('Amplitude');
grid;

subplot(3, 2, 6);
plot(t(startplot:endplot), real(Signal1(startplot:endplot)));
title('Unmodulated Carrier');
xlabel('Time (seconds)');
ylabel('Amplitude');
grid;