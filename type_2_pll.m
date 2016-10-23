%This program performs a sampled time domain simulation of an
%analog phase locked loop (Type 2 with multiplier and low-pass filter)
%
%Written by Aaron Scher
%
%The reference signal is a simple sinusoid. The output should be a sinusoid
%that tracks the frequency of the reference signal after a certain
%start up time.

%Clear work space variables:
clc
clear all

%User inputs:
f0 = 1E6; %Frequency of reference signal [Hz]
phase_ref = 0; %Phase of reference signal [radians]
fVCO = 1.1E6; %free running oscilating frequency of VCO [Hz]
KVCO = 0.5E6; %Gain of VCO (i.e. voltage to frequency transfer coefficient) [Hz/V]
G1 = 0.4467; %Proportional gain term of PI controller
G2 = 1.7783e+05; %Integral gain term of PI controller
fs = 100E6; %Sampling frequency [Hz]
NF = 2000; %Number of samples in simulation
fc = 0.2E6; %Cut-off frequency of low-pass filter (after the multiplier) [Hz]
filter_coefficient_num = 100; %Number of filter coefficeints of low-pass filter

%Start!
b = fir1(filter_coefficient_num,fc/(fs/2)); %design FIR filter coefficients
Ts = 1/fs; %sampling period
t_vec = [0:Ts:(NF-1)*Ts]; %time vector

VCO = zeros(1,NF); %initialize VCO signal array
phi = zeros(1,NF); %initialize VCO angle (phi) array
error = zeros(1,NF); %initialize error array
Int_error = zeros(1,NF); %initialize error array
reference = sin(2*pi*f0*t_vec+phase_ref); %input signal

for n=2:NF
    t=(n-2)*Ts; %Current time (start at t = 0 seconds)

    error_mult(n) = reference(n)*VCO(n-1);%multiply VCO x Signal input to get raw error signal

    %--------Low pass filter the raw error signal:
    for m=1:length(b)
        if n-m+1>=1
           error_array(m) = error_mult(n-m+1);
        else
           error_array(m) = 0;
        end
    end
    error(n) = sum(error_array.*(b));
  %----------------------------------
 
  % Process filtered error signal through PI controller:
           Int_error(n) = Int_error(n-1)+G2*error(n)*Ts;
           PI_error(n) = G1*error(n)+Int_error(n);
  %

    phi(n) = phi(n-1)+2*pi*PI_error(n)*KVCO*Ts; %update the phase of the VCO
    VCO(n) = sin(2*pi*fVCO*t+phi(n)); %compute VCO signal
end

% Plot VCO (output) and reference (input) signals:
figure(1)
plot(t_vec,reference,t_vec,VCO)
title('Plot of input and output signals','FontSize',12)
xlabel('time [s]','FontSize',12)
legend('Input','Output')

%Plot error signal:
figure(2)
plot(t_vec,error)
title('Error signal','FontSize',12)
xlabel('time [s]','FontSize',12) 