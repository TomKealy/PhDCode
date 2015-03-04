function SNR=estimate_SNR(estim_s,true_s)

err = true_s - estim_s;
SNR = 10*log10(sum(abs(true_s).^2)/sum(abs(err).^2));
