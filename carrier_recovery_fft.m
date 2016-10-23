clear all;
close all;

N=10000;
M=20;
Ts = .0001;
time = Ts*(N*M-1);
t = 0:Ts:time;
%m = randi([0 ,1], 1, N);
m = randi([1, 4], 1 ,N)*2 - 3; % pammod(m, 4);
mup = zeros(1, N*M);
mup(1:M:end) = m;
ps = hamming(M);
s = filter(ps, 1, mup);
f0 = 2000;
phoff = -1.0;
c = cos(2*pi*f0*t+phoff);
rsc = s.*c;
rlc = (s+1).*c;

fftrlc = fft(rlc);
fftrsc = fft(rsc);
[m, imax] = max(abs(fftrlc(1:end/2)));
ssf = (0:length(t))/(Ts*length(t));
freqL = ssf(imax);
phaseL = angle(fftrlc(imax));

r = rsc;
q = r.^2;
fl = 500;
ff = [0, .38, .39, .41, .42, 1];
fa = [0 0 1 1 0 0];
h = firpm(fl, ff, fa);
rp = filter(h, 1, q);

fftrBPF = fft(rp);
[mhat, imaxhat] = max(abs(fftrBPF(1:end/2)));
ssx = (0:length(rp))/(Ts*length(rp));
freqS = ssf(imax);
phasep = angle(fftrBPF(imax));
[IR, f] = freqz(h, 1, length(rp), 1/Ts);
[mi, im] = min(abs(f - freqS));
phaseBPF = angle(IR(im));
phaseS = mod(phasep - phaseBPF, pi);