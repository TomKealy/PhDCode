function y = f5252_0(x,SNRb,b)   		
y=Q(-sqrt(2)*x-sqrt(b*SNRb)).^(2^b-1).*exp(-x.*x); % Eq.(5.2.52)
