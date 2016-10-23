%dc040e03.m  for Example 4.3 (Non-uniform Quantization of sin t)
clear, clf
t=[0:0.05:7]; x=sin(t);
b0=-1; bN=1; xmax=max(abs([b0 bN]))
N=16; delta=(bN-b0)/N; % Width of the quantization interval
b=b0+[0:N]*delta; % Boundary vector
c=b(1:N)+delta/2; % Midpoint vector
x_q16=adc(x,b,c); % Uniform quantization
bn=mulaw_inv(b/xmax,50)*xmax;
cn=mulaw_inv(c/xmax,50)*xmax;
xn_q16=adc(x,bn,cn); % Non-uniform quantization
subplot(311), plot(t,x,'k', t,x_q16,'b:', t,xn_q16,'r-')
e=x-x_q16;  en=x-xn_q16; % Absolute error
subplot(312), plot(t,e,'b:', t,en,'r-')
tmp= abs(x)+0.01; re=e./tmp;  ren=en./tmp; % Relative error
subplot(313), plot(t,re,'b:', t,ren,'r-')
