function y=compensate_phase(x,phase)
% To compensate the received OFDM symbol x by phase deviation
y = x*exp(-j*phase);