function phase_estimate=phase_from_pilot(outofFFT,pm)
% To estimate the phase offset based on the received pilot symbol
% outofFFT: OFDM symbol obtained from FFT at the RCVR
% pm      : the sign of pilot signals for each OFDM symbol 
%           determined by the PN sequence 
phase_estimate = angle(outofFFT([8 22 44 58])*[pm;-pm;pm;pm]);