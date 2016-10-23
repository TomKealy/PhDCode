%do_dm.m
% simulates the DM (Delta Modulation) system in Fig.4.9
clear, clf
T=0.01; nmax=200; t=[0:nmax-1]*T;  
b0=-0.1; bN=0.1; N=2;
for itr=1:2
   delta=(bN-b0)/N; % Width of the quantization interval
   xqhn=0; yqn=0; % Initialize the (memory) buffers.
   for n=1:nmax
      x(n)= sin(6*(n-1)*T); % +.2*(rand-0.5);
      xqh(n)=xqhn; % Predictor(delay) output
      d(n)=x(n)-xqh(n); % Comparator input
      dq(n)= sign(d(n))*delta/2; % Comparator
      %with adaptive adjustment of quantization interval width or not 
      if itr==2
        delta=delta*50^(dq(n)*dqn); % Eq.(4.4.1)
      end 
      dqn=dq(n); xqhn=xqh(n)+dqn;  % Store/Delay
      yq(n)=yqn+dq(n);             % Receiver(DM decoder)
      yqn=yq(n);                   % Store/Delay
   end
   subplot(4,1,2*itr-1)
   plot(t,x,'k-'), hold on
   stairs(t,d,'b-'), stairs(t,dq,'r-'), stairs(t,yq,'k-')
   subplot(4,1,2*itr)
   plot(t,x-yq,'r-')
   sum_of_squared_error=(x-yq)*(x-yq)'
end 
