function [taps, taps_weights, taps_norm] = rayleigh_response(no_of_taps, averaging_interval, tap_weights_db, tap_delays):
%no_of_taps = 4;
%averaging_interval = 2000;
%tap_weights_db = [0 -5 -10 -15];
%tap_delays = [0 5e-6 10e-6 15e-6];

% Conversion of dBs scale into linear one
tap_weights_ln = 10.^(tap_weights_db/10);

% Calculation of normalization factor to make the average power of output
% signal to 1
norm_fact = sqrt(sum(tap_weights_ln));

% Generating a Rayleigh Random Variable with average power 1
taps = 1/sqrt(2)*(randn(averaging_interval, no_of_taps) + j*randn(averaging_interval,no_of_taps));

% Making average power of different taps to specified value
taps_weights = taps.*repmat(sqrt(tap_weights_ln), averaging_interval,1);

% Applying to normalization factor
taps_norm = 1/norm_fact*taps_weights; 
end