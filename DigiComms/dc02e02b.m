%dc02e02b.m
N=8000; NB=40; % the number of samples and bins
g_noise=randn(1,N); % an 1xN random vector with N(0,1)
subplot(223), hist(g_noise,NB) % histogram having NB bins
m=1; sigma=1/2; % average and sigma=sqrt(deviation)
g_noise1=g_noise/2+1; % an 1xNB random vector with N(m,sigma^2)
subplot(224), hist(g_noise1,NB) % histogram
[ns,cs]=hist(g_noise1,NB); dx=cs(2)-cs(1); % bin width
x=m+[-500:500]*(sigma/100); % the range on the x-axis
f=exp(-(x-m).^2/(2*sigma^2))/sqrt(2*pi)/sigma; % Eq.(2.1.25)
hold on, plot(x,N*dx*f,'r')
sum(g_noise1)/N-m, sum((g_noise1-m).^2)/N-sigma^2
