%do_square_filter_clock.m - Square-Law Bit Synchronizer
% warning off MATLAB:conversionToLogical
clear, clf
N=100; Ts=1; T=Ts/N; % Symbol time and sampling time
t=(-3*Ts:T:3*Ts)+eps;
rof=0.5; t1=rof*t/Ts;
% Impulse response of a raised-cosine filter
g=sinc(t/Ts).*cos(pi*t1)./(1-4*t1.^2); % Eq.(P2.4.14)
Ng=length(g); ch_input= zeros(1,Ng+N-1); % lch=length(ch_input);
lw=10*N; ws=zeros(4,lw); % signal buffer
tt=[0:lw-1]*T;
% Design of Digital Bandpass Elliptic Filter 
w0=???????; % Center frequency of the BPF to design
wl=w0-0.05; wu=w0^2/wl; % Lower/Upper 3dB frequencies
flfu=2/T*tan([wl wu]*T/2)*T/pi; % Prewarping
Nord=3; Rp=2; As=40; [B,A]=ellip(Nord,Rp,As,flfu);
NA=length(A); NB=length(B);
N_ITER=80; % Number of symbols for simulation
BNRZ_or_URZ=1; test=1; sigma=0.1; 
w0=zeros(max(NA,NB)-1,1); % Zero initial state for the BPF
for n=1:N_ITER
   if BNRZ_or_URZ==0
     a(n)= (rand>0.5)*2-1; % Bipolar NRZ
     ch_input= [ch_input(N+1:end) a(n)*ones(1,N)];
    else 
     a(n)=rand>0.5; % Unipolar RZ
     ch_input=[ch_input(N+1:end) a(n)*[2*ones(1,N/2) zeros(1,N/2)]];
   end
   for m=1:N % Operation per one symbol interval
      r=ch_input(end-N+m-Ng+1:end-N+m)*g'*T +sigma*(rand-0.5);
      ws(1,:)=[ws(1,2:end) r]; % Channel output
      ws(2,:)=[ws(2,2:end) ws(1,end)??]; % Square-law device
      [y,w0]=filter(B,A,ws(2,end),w0); % BPF output
      ws(3,:)=[ws(3,2:end) y(end)];
   end
   ws(4,:) = (ws(3,:)>0);  % Comparator output
   if test>0&&n>N_Iter-20
     for i=1:4, subplot(4,1,i), plot(tt,ws(i,:)), pause;  end
   end
end
