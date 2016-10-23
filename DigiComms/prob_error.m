function p=prob_error(SNRbdB,signaling,b,opt1,opt2)
% Finds the symbol/bit error probability for given SNRbdB (Table 7.1)
%Copyleft: Won Y. Yang, wyyang53@hanmail.net, CAU for academic use only
if nargin<5, opt2 = 'coherent'; end  % opt2='coherent' or 'noncoherent'
if nargin<4, opt1 = 'SER'; end % opt1='SER' or 'BER'
M=2^b; SNRb=10.^(SNRbdB/10);  NSNR=length(SNRb);
if signaling(1:3)=='ASK' % ASK (PAM)
  if lower(opt2(1))=='c' % ASK coherent --> Eq.(7.1.5)
    for i=1:NSNR, p(i)=2*(M-1)/M*Q(sqrt(3*b*SNRb(i)/(M^2-1))); end
    if lower(opt1(1))=='b', p = p/b; end
   else  % ASK non-coherent -->  Eq.(7.1.22)
    if b==1,  for i=1:NSNR, p(i)=exp(-SNRb(i)/4)/2; end;  end   
  end 
elseif signaling(1:3)=='FSK' 
  tmp=M/2/(M-1);
  f5251_=inline('Q(-sqrt(2)*x-sqrt(b*SNRb)).^(2^b-1)','x','SNRb','b');
if lower(opt2(1))=='c' % FSK coherent 
    if b==1
      for i=1:NSNR, p(i)=Q(sqrt(SNRb(i)/2)); end   %Eq.(7.2.9)
     else       
      for i=1:NSNR 
         p(i) = 1-Gauss_Hermite(f5251_,10,SNRb(i),b)/sqrt(pi);
      end   
     end 
  else % FSK non-coherent
    for i=1:NSNR
       p(i)=(M-1)/2*exp(-b*SNRb(i)/4);  tmp1=M-1;
       for m=2:M-1
          tmp1=-tmp1*(M-m)/m; 
          p(i)=p(i)+tmp1/(m+1)*exp(-m*b*SNRb(i)/2/(m+1)); % Eq.(7.2.19)
       end       
    end  
  end 
  if lower(opt1(1))=='b'&b>1, p = p*tmp; end
elseif signaling(1:3)=='PSK'  % Eq.(7.3.7)
  for i=1:NSNR, p(i)=(1+(b>1))*Q(sqrt(b*SNRb(i))*sin(pi/M));  end 
  if lower(opt1(1))=='b'&b>1, p = p/b; end
elseif signaling(1:3)=='DPS' % DPSK
  if b==1  % Eq.(7.4.8)
    for i=1:NSNR, p(i)=2*Q(sqrt(b*SNRb(i)/2)*sin(pi/M)); end 
   else  
    for i=1:NSNR, p(i)=exp(-SNRb(i)/2)/2; end %Eq.(7.4.9)
    if lower(opt1(1))=='b', p = p/b; end
   end
 elseif signaling(1:3)=='QAM'
  L=2^(ceil(b/2)); N = M/L;
  for i=1:NSNR  
     tmpL = 1-2*(L-1)/L*Q(sqrt(3*b/2/(L^2-1)*SNRb(i)));
     tmpN = 1-2*(N-1)/N*Q(sqrt(3*b/2/(N^2-1)*SNRb(i)));
     p(i) = 1-tmpL*tmpN; % Eq.(7.5.6)
  end  
  if lower(opt1(1))=='b'&b>1, p = p/b; end
end
