%Priamry User Signal
path(path,'./Optimization')


% W=6000 bin singal divided into L=500 subbands of bandwidth B=12 bins
% (n.b. all in Fourier domain).

bandwidth_of_spectrum=600;
num_subands=200;
B=bandwidth_of_spectrum/num_subands;
spectrum = zeros(bandwidth_of_spectrum,1);
K = 0.2*bandwidth_of_spectrum; % 2% sparsity

positions = randi(bandwidth_of_spectrum,[1,K]);%generate random spikes for signal

signal = zeros(1,bandwidth_of_spectrum);
signal(positions) = 1;

%signal_per_suband = zeros(500,12);


%plot(signal);

%Secondary User system

%num_nodes = 200;

h = sqrt(1/2)*(randn(bandwidth_of_spectrum,1)+1j*randn(bandwidth_of_spectrum,1));
h = diag(h);

ts = 2*bandwidth_of_spectrum;
fd = 1/(20*ts);
chan = rayleighchan(1e-5,130);

y = filter(chan,signal);