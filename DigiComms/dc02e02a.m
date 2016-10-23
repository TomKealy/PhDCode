%dc02e02a.m
N=8000; NB=40; % the number of samples and bins
u_noise=rand(1,N); % an 1xN random vector with U(0,1)
subplot(221), hist(u_noise,NB) % histogram having NB bins
u_noise1=2*u_noise-1; % an 1xN random vector with U(-1,1)
subplot(222), hist(u_noise1,NB) % histogram
