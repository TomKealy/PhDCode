%test_DSTBC_H4_PSK.m
clear, clf
ITR=10; MaxIter=100; 
sq2= sqrt(2); sq3= sqrt(3); sq6=sqrt(6); 
M=4; phs= 2*pi*[0:M-1]/M; MPSKs= exp(j*phs); % 4-PSK symbols
h=[0.95*exp(j*0.5) exp(-j*0.4) 1.02*exp(-j*0.2) 0.98*exp(j*0.1)]; % Channel 
dh= 0.01*[exp(-j*0.1) -exp(j*0.2) exp(j*0.1) -exp(-j*0.2)]; % Channel variation
Amp_noise = 0.2; % Amplitude of additive Gaussian noise 
ner=0;
for itr1=1:ITR
   S= eye(4); % S0=I;
   r= h*S + Amp_noise*[randn(1,4)+j*randn(1,4)]; % Eq.(P9.10.3)
   x=[]; xh=[]; % Initialize the sequence of transmitted signal
   nn=[-3:0];
   for n=1:MaxIter
      nn= nn+4;
      xn= MPSKs(randint(1,3,M)+1); %  +-1+-j 
      x= [x xn]; % Sequence of transmitted symbols
      x1=xn(1)/sq3; x2=xn(2)/sq3; x3=xn(3)/sq6; x3c=x3';
      x1R=real(x1); x1I=imag(x1); x2R=real(x2); x2I=imag(x2); 
      X=[x1 -x2' x3c  x3c; x2  x1' x3c  -x3c; 
         x3 x3 -x1R+j*x2I x2R+j*x1I; x3 -x3 -x2R+j*x1I -x1R-j*x2I]; %Eq.(P9.10.2)
      S = S*X; % Encoded signal Eq.(P9.10.1)
      r0=r; % Previously received signal
      rs = (h+dh*n)*S; % Received signal Eq.(P9.10.3)
      noise = Amp_noise*[randn(1,4)+j*randn(1,4)]; % Noise components
      r = rs + noise; % Received signal
      r034_=(r0(?)-r0(?))'; r034_r3=r034_*r(3); r034_r4=r034_*r(4);
      r034 = r0(?)+r0(?);   r034r3=r034*r(3)';  r034r4=r034*r(4)';
      xhn1= r0(?)'*r(?)+r0(2)*r(2)'+(-r034_r3+r034_r4-r034r3-r034r4)/2;
      xhn2= r0(2)'*r(?)-r0(?)*r(2)'+(r034_r3+r034_r4-r034r3+r034r4)/2;
      xhn3= (r0(3)'*(r(?)+r(2))+r0(4)'*(r(?)-r(2))+(r0(?)+r0(2))*r(3)'+(r0(?)-r0(2))*r(4)')/sq2;
      xhns= [xhn1 xhn2 xhn3]*(sq3/(r0*r0')); % Decoded signal Eq.(P9.10.4)
      xhn= PSK_slicer(xhns,M);
      xh= [xh xhn]; 
   end
   ner= ner + sum(abs(x-xh)>0.1);
end
SER = ner/(MaxIter*3*ITR) % Symbol error rate