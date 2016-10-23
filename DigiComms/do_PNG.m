%do_PNG.m
% generates m-sequences, Gold sequences, Kasami sequences,
% and then finds their autocorrelations and crosscorrelations
clear
gm=[5 6]; fm=[1 4 5 6]; % Generator polynomials for two LFSRs
m=max(gm); KC=0; x0=[zeros(1,m-1) 1];
g=PNG(gm,KC,x0); gb=g*2-1; % An m-sequence and its bipolarized version
f=PNG(fm,KC,x0); fb=f*2-1; % Another m-sequence and its bipolarized ver
phi_g=corr_circular(gb,gb); % Circular autocorrelation
phi_gf=corr_circular(gb,fb) % Circular crosscorrelation
% 3 values of crosscorrelation by Eq.(10.1.4)
Nc=2^m-1; tm=(1+2^(0.5*(m+1+(mod(m,2)==0)))); 
[unique(sort(phi_gf)); [-tm -1 tm-2]/Nc]
% Two Gold sequences
cG1=rem(g+f,2); cG1b=cG1*2-1; % A Gold sequence 
k=1; % Sequence index in Eg.(10.1.5)
cG2=rem(g+rotate_l(f,k),2); cG2b=cG2*2-1; % Another Gold sequence
phi_G11=corr_circular(cG1b,cG1b); % Circular autocorrelation
phi_G12=corr_circular(cG1b,cG2b); % Circular crosscorrelation
clf, nn=[0:Nc-1]; 
subplot(421), plot(nn,phi_g)
title('Autocorrelation of an m-sequence'), axis([0 Nc -0.3 1])
subplot(422), plot(nn,phi_gf), axis([0 Nc -0.3 1])
title('Crosscorrelation between a preffered pair of m-sequences')
subplot(423), plot(nn,phi_G11), axis([0 Nc -0.3 1])
title('Autocorrelation of a Gold sequence')
subplot(424), plot(nn,phi_G12), axis([0 Nc -0.3 1])
title('Crosscorrelation between two Gold sequences')
% Conversion into Generator polynomial used by MATLAB/Simulink
[gM,m]=gm2gM(gm); [fM,m]=gm2gM(fm); 
sim('PN_generator_sim',Nc)
PN_s=PN(1:Nc).'; [g; PN_s]   % Agree with the Simulink result?
GN_s=GN(1:Nc).'; [cG2; GN_s] % Agree with the Simulink result?
% Prepare for generating Kasami sequences
d=2^(m/2)+1; n=d; g_d=g(n);
while length(g_d)<Nc, n=mod(n-1+d,Nc)+1; g_d=[g_d g(n)]; end
d1=2^(m/2+1)+1; n=d1; g_d1=g(n);
while length(g_d1)<Nc, n=mod(n-1+d1,Nc)+1; g_d1=[g_d1 g(n)]; end
% Three Kasami sequences
k=1; cK1=rem(g+g_d,2); cK2=rem(g+rotate_l(g_d,k),2); cK3=rem(g+g_d1,2);
KN_s1=KN1(1:Nc).'; KN_s2=KN2(1:Nc).'; KN_s3=KN3(1:Nc).';
cK1b=cK1*2-1; cK2b=cK2*2-1; cK3b=cK3*2-1;
cKs1b=KN_s1*2-1; cKs2b=KN_s2*2-1; cKs3b=KN_s3*2-1;
phi_K11=corr_circular(cK1b,cK1b); % Circular autocorrelation
phi_K12=corr_circular(cK1b,cK2b); % Circular crosscorrelation
phi_K13=corr_circular(cK1b,cK3b); % Circular crosscorrelation
[max(abs(phi_K12)) d/Nc max(abs(phi_K13))  d1/Nc]
tn=1+2^(m/2+1); sn=(tn+1)/2; [-tn -sn -1 sn-2 tn-2]/Nc
unique(sort(phi_K12)), unique(sort(phi_K13)) 
subplot(425), plot(nn,phi_K11)
title('Autocorrelation of a Kasami sequence'), axis([0 Nc -0.3 1]), subplot(426), plot(nn,phi_K12)
title('Crosscorrelation between Kasami sequences'), axis([0 Nc -0.3 1])