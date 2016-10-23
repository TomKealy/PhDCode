%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creative Commons
% Attribution-Noncommercial 2.5 India
% You are free:
% to Share — to copy, distribute and transmit the work
% to Remix — to adapt the work
% Under the following conditions:
% Attribution. You must attribute the work in the manner 
% specified by the author or licensor (but not in any way 
% that suggests that they endorse you or your use of the work). 
% Noncommercial. You may not use this work for commercial purposes. 
% For any reuse or distribution, you must make clear to others the 
% license terms of this work. The best way to do this is with a 
% link to this web page.
% Any of the above conditions can be waived if you get permission 
% from the copyright holder.
% Nothing in this license impairs or restricts the author's moral rights.
% http://creativecommons.org/licenses/by-nc/2.5/in/

% Script for simulating binary phase shift keyed transmission and
% reception and compare the simulated and theoretical bit error
% probability
% Checked for proper operation with Octave Version 3.0.0
% Author	: Krishna
% Email		: krishna@dsplog.com
% Version	: 1.0
% Date		: 5 August 2007
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear
N = 10^6 % number of bits or symbols
rand('state',100); % initializing the rand() function
randn('state',200); % initializing the randn() function

% Transmitter
ip = rand(1,N)>0.5; % generating 0,1 with equal probability
s = 2*ip-1; % BPSK modulation 0 -> -1; 1 -> 1 
n = 1/sqrt(2)*[randn(1,N) + j*randn(1,N)]; % white gaussian noise, 0dB variance 
Eb_N0_dB = [-3:10]; % multiple Eb/N0 values

for ii = 1:length(Eb_N0_dB)
   % Noise addition
   10^(-Eb_N0_dB(ii)/20)
   y = s + 10^(-Eb_N0_dB(ii)/20)*n; % additive white gaussian noise

   % receiver - hard decision decoding
   ipHat = real(y)>0;

   % counting the errors
   nErr(ii) = size(find([ip- ipHat]),2);

end

simBer = nErr/N; % simulated ber
theoryBer = 0.5*erfc(sqrt(10.^(Eb_N0_dB/10))); % theoretical ber

% plot
close all
figure
semilogy(Eb_N0_dB,theoryBer,'b.-');
hold on
semilogy(Eb_N0_dB,simBer,'mx-');
axis([-3 10 10^-5 0.5])
grid on
legend('theory', 'simulation');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
title('Bit error probability curve for BPSK modulation');

