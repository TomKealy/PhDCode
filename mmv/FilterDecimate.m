function [Samples] = FilterDecimate(Sig,decfactor,FIR_h);
% low pass filtering and downsampling
f_lpf = fft(FIR_h);
f_sig = fft(Sig);
filter_out = ifft(f_lpf.*f_sig);
Samples = downsample(filter_out,decfactor);
