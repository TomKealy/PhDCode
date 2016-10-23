%test_Rayleigh_fading.m
% To see Rayleigh fading effect (Fig. 2.8)
clear, clf
pi2=2*pi; N=5; % Number of paths
Ns=10000; NB=40; % Numbers of samples and bins
rand('twister',5489) % return rand() to its default initial state
for m=1:Ns
   An=rand(1,N); thn=pi2*rand(N,1);
   rI(m)=An*cos(thn); rQ(m)=An*sin(thn);
end
B=5; rQ = rQ(find(abs(rQ)<=B)); % rQ within the boundary
subplot(311), hist(rQ,NB), hold on
[ns,cs]=hist(rQ,NB); 
dz=cs(2)-cs(1);
z= -B:0.01:B; 
mz=0; sgmz2=1; sgmz=sqrt(sgmz2);
fz= exp(-(z-mz).^2/2/sgmz2)/sqrt(2*pi)/sgmz; % normal (Gaussian) pdf
plot(z,length(rQ)*dz*fz,'r')
As=[0 2];  zmin=0; zmax=5; fzmax=700; 
for i=1:length(As)
A=As(i); % the amplitude of LOS component through a direct path
z= sqrt((A+rI).^2+rQ.^2);
zz=zmin:0.02:zmax; % the range on the z(amplitude)-axis
z=z(find(z>zmin)); z=z(find(z<zmax));
subplot(311+i)
hist(z,NB), hold on
[ns,cs]=hist(z,NB); dz=cs(2)-cs(1);
fz= Rice_pdf(zz,A,1); 
plot(zz,length(z)*dz*fz,'r')
end  
