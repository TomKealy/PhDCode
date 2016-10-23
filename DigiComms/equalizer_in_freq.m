function Y_eq = equalizer_in_freq(Y,H_est)
H_est(find(abs(H_est)<1e-6))=1;  Y_eq = Y./H_est;  % Eq.(11.4.4)