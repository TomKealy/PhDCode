num_samples = 400;

pll_bandwidth = 0.02;
beta = sqrt(pll_bandwidth);

sig_in = zeros(1, num_samples);
phase_in = zeros(1, num_samples);
phase_in(1) = 3.0;
sig_out = zeros(1, num_samples);
phase_error = zeros(1, num_samples);
phase_out = zeros(1, num_samples);
freq_out = zeros(1, num_samples);

for nn = 2:num_samples
  sig_in(nn) = exp(1i*phase_in(nn-1));
  sig_out(nn) = exp(1i*phase_out(nn-1));
  
  phase_error(nn) = arg(sig_in(nn) * conj(sig_out(nn)));
  
  phase_out(nn) =  phase_out(nn-1) + beta*phase_error(nn);
  freq_out(nn) = freq_out(nn-1) + pll_bandwidth*phase_error(nn);
  
  phase_in(nn) = phase_in(nn-1) + freq_in;
  phase_out(nn) = phase_out(nn-1) + freq_out(nn-1);
 end