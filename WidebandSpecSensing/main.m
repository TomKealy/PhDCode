%% Simulation of the following paper 
% A Wideband Spectrum Sensing Method for Cognitive Radio using Sub-Nyquist Sampling
%  IEEE DSP/SPE, Arizona, USA, 2011
clear all;
close all;
clc;
addpath('functions');
%% initial parametrs 
Fs=20;           % Sampling frequency
T=1/Fs;          % Sample time
LL = 1024;        % Length of signal
t = (0:LL-1)*T;  % Time vector
th_norm=20;
%% fft of xci samples
NFFT = 2^nextpow2(LL); % Next power of 2 from length of x 
%% input signal characteristic
%fi=[3.2 7.8 15.6];
fi=[4.8 10.45 15.4];
%fi=[12.5189   16.6941   18.4194];
B=0.9;
F=union(fi-B/2,fi+B/2);
%number of bands
N=length(fi);
ti=[6,13,19]*2;
Ei=[1,1.2,.9,1.3]*4;
Bi=[];
for k=1:2:N*2
    Bi=[Bi F(k+1)-F(k)];
end
B=max(Bi);
%% lebesque measurement
lamda=sum(F(2:2:2*N))-sum(F(1:2:2*N));
%Occupancy 
omega=lamda/Fs;
%% maximom number of active cells
L=floor(Fs/B);
q_max=N+N*ceil(B*L/Fs);
q_min=ceil(N*B*L/Fs);
p=q_max+1;
%% support set of signal
S=[];
FF=ceil(F*L*T);
for h=1:2:N*2
    S=[S, max(FF(h)):FF(h+1)];
end    
S=unique(S)-1;
%% pick a sample pattern with SFS
C=SFS_C(L,p,S+1,Fs);
%% A matrix
k=1:L;
A=1/(L*T)*exp(1i*2*pi*C'*(k-1)/L);
As=A(:,S+1);
%% check for low condition number 
while cond(As)> 10
    C=pickn(L,p);
    k=1:L;
    A=1/(L*T)*exp(1i*2*pi*C'*(k-1)/L);
    As=A(:,S+1);    
end
%% input signal generation 
x=zeros(1,LL);
for n=1:N
    x=x+sqrt(Ei(n)*Bi(n)) * sinc(Bi(n)*(t-ti(n))).* exp(1i*2*pi*fi(n)*t);
end
sigma=1/256;
w=sigma*(randn(1,LL)+1i*randn(1,LL))/sqrt(2);
x=x+w;
%% fft of original signal
Xf=fft(x,NFFT);
%% multi coset sampling
xci=zeros(p,LL);
for k=1:p
    for m=C(k)+1:L:LL-C(k)
        xci(k,m)=x(m);
    end
end
%% fir filter
Ntap=191*2+1;
hr=fircls1(Ntap,1/L,0.02,0.008);
n=0:Ntap;
h=hr.*exp(i*pi*n/L);
xci_h=zeros(p,LL);
%% filter coset samples
d1=(Ntap+1)/2;
for k=1:p
    xfilter=conv(xci(k,:),h);
    xci_h(k,:)=xfilter(d1+1:end-d1+1);
end
% sigma=1/256;
% w=sigma*(randn(p,LL)+i*randn(p,LL))/sqrt(2);
% xci_h=xci_h+w;
%% R=<x,x> 
R= xci_h * xci_h';
%% eigen value decomposition
% find eigenvalues and eigenvectors
[Us,Gama]=eig(R);
ev=eig(R);
figure;
plot(sort(ev,'descend'),'--s','LineWidth',3,'MarkerEdgeColor','k',...
                                 'MarkerFaceColor','g',...
                                 'MarkerSize',10);
title('Eigenvalues');
xlabel('i','fontsize',12,'fontweight','b');
ylabel('\lambda _i','fontsize',12,'fontweight','b');

box off;
%set(gca,'YTickLabel',[]);
%% find number of active slots with MDL
q_hat=mdl_function(R,LL);
%q_hat=aic_function(R,LL);
%% find location of active slots
%create Un:matrix of eigenvectors correspond to zero eigenvalues
Un=Us(:,1:p-q_hat);
UnA=Un'*A;
for k=1:L
    nUnA(k)=norm(A(:,k)'*Un);
%    nUnA(k)=norm(UnA(:,k),2);
end
figure;
stem(0:L-1,1./nUnA,'--','LineWidth',3,'MarkerEdgeColor','k',...
                                 'MarkerFaceColor','g',...
                                 'MarkerSize',10);
xlabel('k, spectral index','fontsize',12,'fontweight','b');
ylabel('P_M_U(k)','fontsize',12,'fontweight','b');
box off;
axis([0 L-1 0 max(1./nUnA)+10]);
% find q minimom norm of Un A to find spectral support
pmu=1./nUnA;
[pmus,loc]=sort(pmu);
Sr=loc(end-q_hat+1:end);
Sr=sort(Sr)-1;
q=length(Sr);
%title(['Sr=[',num2str(Sr),']','  q_I_C=',num2str(q_hat)]);
%% save only Sr culumns
As=A(:,Sr+1);
pAs=pinv(As);
K_As=cond(As)
%% reconstruction in time domain
t = (0:NFFT-1)*T;  % Time vector
xt=zeros(1,NFFT);
for n=1:NFFT
for m=1:q
    for l=1:p
        xt(n)=xt(n)+pAs(m,l)*xci_h(l,n)*exp(i*2*pi*Sr(m)*(n-1)/L)/T;
    end    
end
end
%% DFT of reconstructed signal
Xtf=fft(xt,NFFT);
%% plots
figure;
subplot(2,1,1);
f = Fs*linspace(0,1,NFFT);
plot(f,abs(Xf),f,abs(Xtf),'r');
axis([0 20 0 abs(max(Xtf))]);
yb=1:80;
for gr=Fs/L:Fs/L:Fs
    hold on;
    plot(gr,yb,'g');
    text(gr-2*Fs/(3*L),yb(end)-10,num2str(gr*L/Fs));
end    
xlabel('f');
ylabel('X(f)');
title(['S=[',num2str(S),']']);
%% MSE
MSE=100*norm(abs(x)-abs(xt))/norm(abs(x))
subplot(2,1,2);
plot(Fs*t(1:LL),abs(x(1,1:LL)),Fs*t(1:LL),abs(xt(1,1:LL)),'r');
xlabel('t');
title(['L=',num2str(L),'     p=',num2str(p),'     MSE=',num2str(MSE),'%']);
%%
figure;
subplot(2,1,1)
stem(0:L-1,1./nUnA,'--','LineWidth',3,'MarkerEdgeColor','k',...
                                 'MarkerFaceColor','g',...
                                 'MarkerSize',10);
xlabel('k, spectral index','fontsize',12,'fontweight','b');
ylabel('P_M_U(k)','fontsize',12,'fontweight','b');
box off;
axis([0 L-1 0 max(1./nUnA)+10]);

subplot(2,1,2)
plot(f,abs(Xf),'LineWidth',3);
axis([0,Fs,0,50]);
box off
xlabel('frequency','fontsize',12,'fontweight','b');
ylabel('|X(f)|','fontsize',12,'fontweight','b');
%set(gca,'YTickLabel',[]);
%set(gca,'XTickLabel',[]);
yb=1:80;
for gr=Fs/L:Fs/L:Fs
    hold on;
    plot(gr,yb,'g');
    text(gr-2*Fs/(3*L),yb(end)-10,num2str(gr*L/Fs));
end    
%%
Xf=fft(x,NFFT);
figure;
subplot(2,2,1);
plot(t,abs(x));
axis([0,t(end),0,abs(max(x))]);
box off
xlabel('time');
ylabel('|x(t)|');

subplot(2,2,3);
plot(t,abs(xt),'r');
axis([0,t(end),0,abs(max(xt))]);
box off
xlabel('time');
ylabel('|x_r(t)|');

subplot(2,2,2);
plot(f,abs(Xf));
axis([0,Fs,0,50]);
box off
xlabel('frequency');
ylabel('|X(f)|');


subplot(2,2,4);
plot(f,abs(Xtf),'r');
axis([0,Fs,0,50]);
box off
xlabel('frequency');
ylabel('|X_r(f)|');
