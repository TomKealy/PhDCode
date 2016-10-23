%do_mmse.m to simulate MMSE equalizer
clear
T=0.5; t=-3*T:T/2:3*T;
g=1./(1+(2*t/T).^2); %g(1:6)=g(1:6)/5;
g=[0.027 0.0385 0.0588 0.1  0.2  0.5  1  0.5  0.2  0.1 0.0588 0.0385 0.027];
lg=length(g);  
[gmax,n01]=max(g); n0=n01-1; K=5;
d=n0+floor(K/2); %expected delay of the channel and equalizer by (6.2-12)
Na=1000; a=(rand(1,Na)>0.5)*2-1; %data sequence
sigma=0.2;
y=filter(g,1,a)+sigma*(rand(1,Na)-0.5);
c_MMSEE=mmsee(y,a,K,d) %mmse equalizer coefficients
gn=g+sigma*(rand(size(g))-0.5); % channel response with noise
g_MMSEE =filter(c_MMSEE,1,gn); %mmse output to noisy channel response

figure(1), clf
subplot(611), stem([0:lg-1],g, 'Markersize',5)
title('Original channel impulse response')
set(gca,'fontsize',9), axis([-1 lg -0.3 1.2])
hold on, plot([-1 lg], [0 0])

for d=0:4
  c=mmsee(y,a,K,n0+d)
  ge =filter(c,1,g);  
  subplot(612+d), stem([0:lg-1],ge)
  hold on, plot([-1 lg], [0 0])
  set(gca,'fontsize',9), axis([-1 lg -0.3 1.2])
end

%To compare MMSEE with Zero-forcing Equalizer
%gc=g([n01-(K-1):n01+(K-1)]);
%gc=gc+sigma*(rand(size(gc))-0.5);
d=floor(K/2); c_ZFE=zfe(gn,K,d) % ZFE equalizer coefficients
g_ZFE =filter(c_ZFE,1,gn); % ZFE output to (channel response+noise)

figure(2), clf
subplot(511), stem([0:lg-1],gn, 'Markersize',5)
title('AWGN channel impulse response')
set(gca,'fontsize',9), axis([-1 lg -0.3 1.2])
hold on, plot([-1 lg], [0 0])
subplot(512), stem([0:lg-1],g_ZFE, 'Markersize',5)
title('ZFE response')
set(gca,'fontsize',9), axis([-1 lg -0.3 1.2])
hold on, plot([-1 lg], [0 0])
subplot(513), stem([0:lg-1],g_MMSEE, 'Markersize',5)
title('MMSEE response')
set(gca,'fontsize',9), axis([-1 lg -0.3 1.2])
hold on, plot([-1 lg], [0 0])

%adaptive equalizer
c_ADE= zeros(1,K);  
ch_input= zeros(1,lg);
d=floor(K/2); delta=0.05;
for m=1:Na
   a(m)= (rand>0.5)*2-1;
   ch_input= [a(m) ch_input(1:lg-1)]; 
   noise(m)= sigma*(rand-0.5);
   y(m)= ch_input*g' + noise(m);
   if m>max(K-1,n0+d)
     [c_ADE,se(m)]=ade(c_ADE,y(m:-1:m-K+1),a(m-(n0+d)),delta);
   end
end
g_ADE =filter(c_ADE,1,gn); %ade output to (channel response+noise)
subplot(514), stem([0:lg-1],g_ADE, 'Markersize',5)
title('ADE response')
set(gca,'fontsize',9), axis([-1 lg -0.3 1.2])
hold on, plot([-1 lg], [0 0])
subplot(515), plot(se)
title('Squared error between the ADE output and the desired output (delayed channel input)')
set(gca,'fontsize',9), axis([0 Na 0 2.5])
figure(3), subplot(211), plot(10*log10(se)), subplot(212), semilogy(se,'r')
